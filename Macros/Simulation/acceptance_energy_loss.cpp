#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TApplication.h>
#include <iostream>
#include <array>
#include <cmath>
#include <string>
#include <exception>

#include "catima/catima.h"

using namespace std;

// --- Geometría nominal (valores base) ---
const double detector_distance_nominal = 5.0; // cm
const double detector_width = 20.0;           // cm
const double detector_height = 20.0;          // cm
const double target_radius = 4.0;             // cm
const double target_thickness_cm = 0.41e-4;   // cm
const double backing_thickness_cm = 2.0e-4;   // cm
const double gas_path_nominal_cm = 2.5;       // cm
const double gas_path_cathode_cm = 0.32;      // cm

// ---- Material parameters ----
const double target_uo2_areal_gcm2 = 0.3e-3;
const double target_ceo2_areal_gcm2 = 1.2e-3;
const double backing_areal_gcm2 = 2.5e-4 * 2.70;
const double gas_density_gcm3 = 0.000032;

// --- Configuración geométrica a variar ---
struct GeoConfig {
    double rot_angle_deg;   // ángulo de rotación primed->beam
    double det_dist_cm;     // distancia entre detectores
    double dx_front;        // offset en x del detector frontal
    double dy_front;        // offset en y del detector frontal
    double dx_back;         // offset en x del detector trasero
    double dy_back;         // offset en y del detector trasero
    std::string tag;        // sufijo para nombre de archivo
};

// --- Función TOF a través de una capa ---
double tof_through_layer(double E_in_mevu, double E_loss, int A, double path_cm){
    double E_in = E_in_mevu * A;
    double E_out = E_in - E_loss;

    if (E_in <= 0.0) return 1e9;
    double E_avg = 0.5*(E_in + E_out);
    const double u = 931.494;
    const double c_cm_per_ns = 29.9792458;
    double mass = A * u;
    double speed = sqrt(2.0 * E_avg / mass);
    double v = c_cm_per_ns * speed;
    if (v <= 0.0) return 1e9;
    return path_cm / v;
}

// intersection between the fission fragment trajectory and the plane of the detector
bool line_plane_intersection(const array<double,3>& origin, const array<double,3>& dir,
                             double z_plane, array<double,3>& out){
    double dz = dir[2];
    if (fabs(dz) < 1e-12) return false;
    double t = (z_plane - origin[2]) / dz;
    if (t <= 0.0) return false;
    out[0] = origin[0] + t*dir[0];
    out[1] = origin[1] + t*dir[1];
    out[2] = origin[2] + t*dir[2];
    return true;
}

// check if the fragment hits the detector area
bool in_detector(const array<double,3>& p, bool front,
                 const GeoConfig &cfg){
    double x_high_limit, x_low_limit;
    double dx = front ? cfg.dx_front : cfg.dx_back;
    double dy = front ? cfg.dy_front : cfg.dy_back;

    double x = p[0] - dx;
    double y = p[1] - dy;

    if (front) {
        x_high_limit = detector_width/2.0 + 2.5;
        x_low_limit  = -detector_width/2.0 + 2.5;
    } else {
        x_high_limit = detector_width/2.0 - 2.5;
        x_low_limit  = -detector_width/2.0 - 2.5;
    }
    return (x >= x_low_limit) && (x <= x_high_limit) &&
           (fabs(y) <= detector_height/2.0);
}

// function to rotate from the beam to the detector frame
void primed_to_beam_coordinates(const array<double,3>& p,
                                array<double,3>& out,
                                double rot_angle_deg){
    double cA = cos(TMath::DegToRad()*rot_angle_deg);
    double sA = sin(TMath::DegToRad()*rot_angle_deg);
    double x_b = p[0] * cA + p[2] * sA;
    double y_b = p[1];
    double z_b = -p[0] * sA + p[2] * cA;
    out[0]=x_b; out[1]=y_b; out[2]=z_b;
}

// run different configurations
void run_simulation(const GeoConfig &cfg,
                    const char* infile_base){
    // archivo de salida con sufijo
    std::string outfile = std::string("fission_events_output_uranium_toy_") + cfg.tag + ".root";

    TRandom3 rng(42);
    TFile *fin = TFile::Open(infile_base, "READ");
    if (!fin || fin->IsZombie()){
        cerr << "Cannot open input file: " << infile_base << "\n";
        return;
    }

    TTree *tin = (TTree*) fin->Get("FissionTree");
    if (!tin){
        cerr << "FissionTree not found.\n";
        fin->Close();
        return;
    }

    // Input branches
    Int_t Z1=0, Z2=0;
    Int_t A1post=0, A2post=0;
    Double_t tke_tot=0, tke1=0, tke2=0;
    Double_t neutron_energy=0;

    tin->SetBranchAddress("Z1", &Z1);
    tin->SetBranchAddress("Z2", &Z2);
    tin->SetBranchAddress("A1post", &A1post);
    tin->SetBranchAddress("A2post", &A2post);
    tin->SetBranchAddress("tke", &tke_tot);
    tin->SetBranchAddress("tke1", &tke1);
    tin->SetBranchAddress("tke2", &tke2);

    // Output file & tree
    TFile *fout = TFile::Open(outfile.c_str(), "RECREATE");
    TTree *tout = new TTree("fission_process", "Fission with CATima energy losses and systematics");

    // Output branches 
    Float_t initial_energy_forward_mevu=0, initial_energy_backward_mevu=0;
    Float_t loss_target_forward=0, loss_target_backward=0;
    Float_t loss_backing_forward=0;
    Float_t loss_gas_forward=0, loss_gas_backward=0;
    Float_t final_energy_forward_mevu=0, final_energy_backward_mevu=0;
    Float_t tof_diff=0;
    Float_t theta_deg=0, phi_deg=0, cos_theta=0;
    Float_t x_back_local=0, y_back_local=0, x_front_local=0, y_front_local=0;
    Float_t x_back_world=0, y_back_world=0, x_front_world=0, y_front_world=0;
    Int_t in_Z1=0, in_Z2=0, in_A1=0, in_A2=0;
    Float_t in_tke1=0, in_tke2=0;
    Int_t forward_Z=0, forward_A=0, backward_Z=0, backward_A=0;
    Float_t loss_mylar_cathode_x_f=0, loss_mylar_anode_f=0, loss_mylar_cathode_y_f=0;
    Float_t loss_gas_cathode_x_f=0, loss_gas_cathode_y_f=0;
    Float_t loss_mylar_cathode_x_b=0, loss_mylar_anode_b=0, loss_mylar_cathode_y_b=0;
    Float_t loss_gas_cathode_x_b=0, loss_gas_cathode_y_b=0;
    Float_t tof_backward=0, tof_forward=0, tof_backing=0;

    // configuration of the errors
    Float_t br_rot_angle_deg = cfg.rot_angle_deg;
    Float_t br_det_dist_cm   = cfg.det_dist_cm;
    Float_t br_dx_front      = cfg.dx_front;
    Float_t br_dy_front      = cfg.dy_front;
    Float_t br_dx_back       = cfg.dx_back;
    Float_t br_dy_back       = cfg.dy_back;

    // Set up output branches
    tout->Branch("initial_energy_forward_mevu",&initial_energy_forward_mevu);
    tout->Branch("initial_energy_backward_mevu",&initial_energy_backward_mevu);
    tout->Branch("loss_target_forward",&loss_target_forward);
    tout->Branch("loss_target_backward",&loss_target_backward);
    tout->Branch("loss_backing_forward",&loss_backing_forward);
    tout->Branch("loss_gas_forward",&loss_gas_forward);
    tout->Branch("loss_gas_backward",&loss_gas_backward);
    tout->Branch("final_energy_forward_mevu",&final_energy_forward_mevu);
    tout->Branch("final_energy_backward_mevu",&final_energy_backward_mevu);
    tout->Branch("loss_mylar_cathode_x_f",&loss_mylar_cathode_x_f);
    tout->Branch("loss_mylar_anode_f",&loss_mylar_anode_f);
    tout->Branch("loss_mylar_cathode_y_f",&loss_mylar_cathode_y_f);
    tout->Branch("loss_gas_cathode_x_f",&loss_gas_cathode_x_f);
    tout->Branch("loss_gas_cathode_y_f",&loss_gas_cathode_y_f);
    tout->Branch("loss_mylar_cathode_x_b",&loss_mylar_cathode_x_b);
    tout->Branch("loss_mylar_anode_b",&loss_mylar_anode_b);
    tout->Branch("loss_mylar_cathode_y_b",&loss_mylar_cathode_y_b);
    tout->Branch("loss_gas_cathode_x_b",&loss_gas_cathode_x_b);
    tout->Branch("loss_gas_cathode_y_b",&loss_gas_cathode_y_b);
    tout->Branch("tof_diff",&tof_diff);
    tout->Branch("tof_backing", &tof_backing);
    tout->Branch("theta_deg",&theta_deg);
    tout->Branch("phi_deg",&phi_deg);
    tout->Branch("cos_theta",&cos_theta);
    tout->Branch("x_back_local",&x_back_local);
    tout->Branch("y_back_local",&y_back_local);
    tout->Branch("x_front_local",&x_front_local);
    tout->Branch("y_front_local",&y_front_local);
    tout->Branch("x_back_world",&x_back_world);
    tout->Branch("y_back_world",&y_back_world);
    tout->Branch("x_front_world",&x_front_world);
    tout->Branch("y_front_world",&y_front_world);
    tout->Branch("tof_backward",&tof_backward);
    tout->Branch("tof_forward",&tof_forward);
    tout->Branch("in_Z1",&in_Z1);
    tout->Branch("in_Z2",&in_Z2);
    tout->Branch("in_A1",&in_A1);
    tout->Branch("in_A2",&in_A2);
    tout->Branch("in_tke1",&in_tke1);
    tout->Branch("in_tke2",&in_tke2);
    tout->Branch("forward_Z",&forward_Z);
    tout->Branch("forward_A",&forward_A);
    tout->Branch("backward_Z",&backward_Z);
    tout->Branch("backward_A",&backward_A);

    tout->Branch("rot_angle_deg",&br_rot_angle_deg);
    tout->Branch("det_dist_cm",&br_det_dist_cm);
    tout->Branch("dx_front",&br_dx_front);
    tout->Branch("dy_front",&br_dy_front);
    tout->Branch("dx_back",&br_dx_back);
    tout->Branch("dy_back",&br_dy_back);

    Long64_t nEntries = tin->GetEntries();
    cout << "Entries: " << nEntries << " for config " << cfg.tag << "\n";
    Long64_t accepted = 0, stopped = 0, missed = 0;

    // Materials
    catima::Material uo2_ref1;
    uo2_ref1.add_element(238.03, 92, 1.0);
    uo2_ref1.add_element(16.00, 8, 2.0);
    uo2_ref1.density(7.29);

    catima::Material uo2_ref2;
    uo2_ref2.add_element(238.03, 92, 1.0);
    uo2_ref2.add_element(16.00, 8, 2.0);
    uo2_ref2.density(7.29);

    catima::Material al_ref;
    al_ref.add_element(27.0, 13, 1.0);
    al_ref.density(2.70);

    catima::Material c3f8_ref;
    c3f8_ref.add_element(12.0, 6, 3.0);
    c3f8_ref.add_element(19.0, 9, 8.0);
    c3f8_ref.density(0.000030855);

    catima::Material mylar_cathode_x;
    mylar_cathode_x = catima::get_material(214);
    mylar_cathode_x.thickness(0.00020955);

    catima::Material mylar_anode;
    mylar_anode = catima::get_material(214);
    mylar_anode.thickness(0.00020955);

    catima::Material mylar_cathode_y;
    mylar_cathode_y = catima::get_material(214);
    mylar_cathode_y.thickness(0.00020955);

    catima::Material c3f8_ref_cathode;
    c3f8_ref_cathode.add_element(12.0, 6, 3.0);
    c3f8_ref_cathode.add_element(19.0, 9, 8.0);
    c3f8_ref_cathode.density(0.000030855);

    //detector z planes
    double z_det_back  =  cfg.det_dist_cm / 2.0;
    double z_det_front = -cfg.det_dist_cm / 2.0;

    for (Long64_t i=0; i<nEntries; ++i){
        tin->GetEntry(i);

        in_Z1 = Z1; in_Z2 = Z2;
        in_A1 = static_cast<int>(round(A1post));
        in_A2 = static_cast<int>(round(A2post));
        in_tke1 = static_cast<float>(tke1);
        in_tke2 = static_cast<float>(tke2);

        bool frag1_forward = rng.Uniform() > 0.5;

        forward_Z  = frag1_forward ? Z1 : Z2;
        forward_A  = frag1_forward ? A1post : A2post;
        backward_Z = frag1_forward ? Z2 : Z1;
        backward_A = frag1_forward ? A2post : A1post;
        double forward_tke  = frag1_forward ? tke1 : tke2;
        double backward_tke = frag1_forward ? tke2 : tke1;

        double E_forward_mevu = forward_tke / forward_A;
        double E_backward_mevu = backward_tke / backward_A;

        double phi_r = rng.Uniform(0, 2*M_PI);
        double r = target_radius * (rng.Uniform(0.0, 1.0));
        double tx = r * cos(phi_r);
        double ty = r * sin(phi_r);
        double tz = rng.Uniform(-target_thickness_cm/2.0, target_thickness_cm/2.0);

        array<double,3> origin = {tx, ty, tz};

        double phi = rng.Uniform(0, 2*M_PI);
        double u = rng.Uniform(-1.0, 0.0);
        double theta = acos(u);
        double cos_theta_val = u;
        double sin_theta_val = sin(theta);

        array<double,3> dir_forward = {sin_theta_val*cos(phi), sin_theta_val*sin(phi), cos_theta_val};
        array<double,3> dir_backward = {-dir_forward[0], -dir_forward[1], -dir_forward[2]};

        array<double,3> hit_front, hit_back;
        bool ok_front = line_plane_intersection(origin, dir_forward, z_det_front, hit_front);
        bool ok_back = line_plane_intersection(origin, dir_backward, z_det_back, hit_back);

        if (!ok_front || !ok_back){ missed++; continue; }
        if (!in_detector(hit_front, true, cfg) || !in_detector(hit_back, false, cfg)){
            missed++;
            continue;
        }

        double path_target_backward = (target_thickness_cm/2.0 - origin[2]) / fabs(cos_theta_val);
        double path_target_forward = (target_thickness_cm/2.0 + origin[2]) / fabs(cos_theta_val);

        double uo2_density = target_uo2_areal_gcm2 / target_thickness_cm;
        double target_areal_forward = uo2_density * path_target_forward;
        double target_areal_backward = uo2_density * path_target_backward;
        uo2_ref1.thickness(target_areal_forward);
        uo2_ref2.thickness(target_areal_backward);

        double path_backing_forward = backing_thickness_cm / fabs(cos_theta_val);
        double backing_areal_forward = backing_areal_gcm2 * (path_backing_forward / backing_thickness_cm);
        al_ref.thickness(backing_areal_forward);

        double path_gas_forward = gas_path_nominal_cm / fabs(cos_theta_val);
        double path_gas_backward = gas_path_nominal_cm / fabs(cos_theta_val);
        double gas_areal_forward = gas_density_gcm3 * path_gas_forward;
        c3f8_ref.thickness(gas_areal_forward);
        double gas_areal_backward = gas_density_gcm3 * path_gas_backward;

        double path_gas_cathode_x = gas_path_cathode_cm / fabs(cos_theta_val);
        double gas_areal_cathode_x = gas_density_gcm3 * path_gas_cathode_x;
        c3f8_ref_cathode.thickness(gas_areal_cathode_x);

        double path_mylar = 1.5e-4/fabs(cos_theta_val);
        double mylar_density = 1.39;
        double mylar_areal = mylar_density * path_mylar;
        mylar_cathode_x.thickness(mylar_areal);
        mylar_anode.thickness(mylar_areal);
        mylar_cathode_y.thickness(mylar_areal);

        catima::Projectile p_forward(forward_A, forward_Z);
        p_forward.T = E_forward_mevu;

        catima::Projectile p_backward(backward_A, backward_Z);
        p_backward.T = E_backward_mevu;

        double E_forward_final_mevu = E_forward_mevu;
        double E_backward_final_mevu = E_backward_mevu;
        double loss_t_f = 0, loss_bk_f = 0, loss_g_f = 0;
        double loss_t_b = 0, loss_g_b = 0;
        loss_mylar_cathode_x_f = 0; loss_mylar_anode_f = 0; loss_mylar_cathode_y_f = 0;
        loss_gas_cathode_x_f = 0; loss_gas_cathode_y_f = 0;
        loss_mylar_cathode_x_b = 0; loss_mylar_anode_b = 0; loss_mylar_cathode_y_b = 0;
        loss_gas_cathode_x_b = 0; loss_gas_cathode_y_b = 0;

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_uo2 = catima::calculate(p_forward, uo2_ref1);
            double energy_loss_uo2 = res_uo2.Eloss;
            E_forward_final_mevu = res_uo2.Eout;
            loss_t_f += energy_loss_uo2;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_al = catima::calculate(p_forward, al_ref);
            double energy_loss_al = res_al.Eloss;
            E_forward_final_mevu = res_al.Eout;
            loss_bk_f += energy_loss_al;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_gas = catima::calculate(p_forward, c3f8_ref);
            double energy_loss_gas = res_gas.Eloss;
            E_forward_final_mevu = res_gas.Eout;
            loss_g_f += energy_loss_gas;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_mylar1 = catima::calculate(p_forward, mylar_cathode_x);
            double energy_loss_mylar1 = res_mylar1.Eloss;
            E_forward_final_mevu = res_mylar1.Eout;
            loss_mylar_cathode_x_f += energy_loss_mylar1;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_gas2 = catima::calculate(p_forward, c3f8_ref_cathode);
            double energy_loss_gas2 = res_gas2.Eloss;
            E_forward_final_mevu = res_gas2.Eout;
            loss_gas_cathode_x_f += energy_loss_gas2;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_mylar2 = catima::calculate(p_forward, mylar_anode);
            double energy_loss_mylar2 = res_mylar2.Eloss;
            E_forward_final_mevu = res_mylar2.Eout;
            loss_mylar_anode_f += energy_loss_mylar2;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_gas3 = catima::calculate(p_forward, c3f8_ref_cathode);
            double energy_loss_gas3 = res_gas3.Eloss;
            E_forward_final_mevu = res_gas3.Eout;
            loss_gas_cathode_y_f += energy_loss_gas3;
        }

        if (E_forward_final_mevu > 0.001) {
            p_forward.T = E_forward_final_mevu;
            auto res_mylar3 = catima::calculate(p_forward, mylar_cathode_y);
            double energy_loss_mylar3 = res_mylar3.Eloss;
            E_forward_final_mevu = res_mylar3.Eout;
            loss_mylar_cathode_y_f += energy_loss_mylar3;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_uo2_b = catima::calculate(p_backward, uo2_ref2);
            double energy_loss_uo2_b = res_uo2_b.Eloss;
            E_backward_final_mevu = res_uo2_b.Eout;
            loss_t_b += energy_loss_uo2_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_gas_b = catima::calculate(p_backward, c3f8_ref);
            double energy_loss_gas_b= res_gas_b.Eloss;
            E_backward_final_mevu = res_gas_b.Eout;
            loss_g_b += energy_loss_gas_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_mylar1_b = catima::calculate(p_backward, mylar_cathode_x);
            double energy_loss_mylar1_b = res_mylar1_b.Eloss;
            E_backward_final_mevu = res_mylar1_b.Eout;
            loss_mylar_cathode_x_b += energy_loss_mylar1_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_gas2_b = catima::calculate(p_backward, c3f8_ref_cathode);
            double energy_loss_gas2_b = res_gas2_b.Eloss;
            E_backward_final_mevu = res_gas2_b.Eout;
            loss_gas_cathode_x_b += energy_loss_gas2_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_mylar2_b = catima::calculate(p_backward, mylar_anode);
            double energy_loss_mylar2_b = res_mylar2_b.Eloss;
            E_backward_final_mevu = res_mylar2_b.Eout;
            loss_mylar_anode_b += energy_loss_mylar2_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_gas3_b = catima::calculate(p_backward, c3f8_ref_cathode);
            double energy_loss_gas3_b = res_gas3_b.Eloss;
            E_backward_final_mevu = res_gas3_b.Eout;
            loss_gas_cathode_y_b += energy_loss_gas3_b;
        }

        if (E_backward_final_mevu > 0.001) {
            p_backward.T = E_backward_final_mevu;
            auto res_mylar3_b = catima::calculate(p_backward, mylar_cathode_y);
            double energy_loss_mylar3_b = res_mylar3_b.Eloss;
            E_backward_final_mevu = res_mylar3_b.Eout;
            loss_mylar_cathode_y_b += energy_loss_mylar3_b;
        }

        if (E_forward_final_mevu <= 0.001 || E_backward_final_mevu <= 0.001) {
            stopped++;
            continue;
        }

        initial_energy_forward_mevu = E_forward_mevu;
        initial_energy_backward_mevu = E_backward_mevu;
        loss_target_forward = loss_t_f;
        loss_target_backward = loss_t_b;
        loss_backing_forward = loss_bk_f;
        loss_gas_forward = loss_g_f;
        loss_gas_backward = loss_g_b;
        final_energy_forward_mevu = E_forward_final_mevu;
        final_energy_backward_mevu = E_backward_final_mevu;

        double tof_target_forward = tof_through_layer(E_forward_mevu, loss_t_f, forward_A, path_target_forward);
        double tof_backing_forward = tof_through_layer(E_forward_mevu - loss_t_f/forward_A, loss_bk_f, forward_A, path_backing_forward);
        tof_backing = tof_backing_forward;
        double tof_gas_forward = tof_through_layer(E_forward_mevu - loss_t_f/forward_A - loss_bk_f/forward_A, loss_g_f, forward_A, path_gas_forward);
        double tof_mylar_forward = tof_through_layer(E_forward_mevu - loss_t_f/forward_A - loss_bk_f/forward_A - loss_g_f/forward_A, loss_mylar_cathode_x_f, forward_A, path_mylar);
        tof_forward = tof_target_forward + tof_backing_forward + tof_gas_forward + tof_mylar_forward;

        double tof_target_backward = tof_through_layer(E_backward_mevu, loss_t_b, backward_A, path_target_backward);
        double tof_gas_backward = tof_through_layer(E_backward_mevu - loss_t_b/backward_A, loss_g_b, backward_A, path_gas_backward);
        double tof_mylar_backward = tof_through_layer(E_backward_mevu - loss_t_b/backward_A - loss_g_b/backward_A, loss_mylar_cathode_x_b, backward_A, path_mylar);
        tof_backward = tof_target_backward + tof_gas_backward + tof_mylar_backward;
        tof_diff = tof_forward - tof_backward;

        x_front_local = hit_front[0]; y_front_local = hit_front[1];
        x_back_local = hit_back[0];   y_back_local = hit_back[1];

        array<double,3> hf_beam, hb_beam;
        primed_to_beam_coordinates(hit_front, hf_beam, cfg.rot_angle_deg);
        primed_to_beam_coordinates(hit_back,  hb_beam, cfg.rot_angle_deg);
        x_front_world = hf_beam[0]; y_front_world = hf_beam[1];
        x_back_world  = hb_beam[0]; y_back_world  = hb_beam[1];

        array<double,3> dir_forward_beam;
        primed_to_beam_coordinates(dir_forward, dir_forward_beam, cfg.rot_angle_deg);

        double dx = dir_forward_beam[0];
        double dy = dir_forward_beam[1];
        double dz = dir_forward_beam[2];
        double norm = sqrt(dx*dx + dy*dy + dz*dz);
        if (norm > 0) {
            double ct = dz/norm;
            cos_theta = ct;
            theta_deg = acos(ct) * 180.0 / M_PI;
            phi_deg = atan2(dy, -dx) * 180.0 / M_PI;
        }

        accepted++;

        tout->Fill();
        if (i%1000 == 0){
            std::cout << " accepted: " << accepted
              << ", stopped: " << stopped
              << ", missed: " << missed << "\n";

        }
    }

    tout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Config " << cfg.tag
              << " accepted: " << accepted
              << ", stopped: " << stopped
              << ", missed: " << missed << "\n";
}
int main(){
    const char* infile = "/Users/nico/root-on-vscode/fission_events.root";

    // variaciones de ±5 % en pasos de 1 % para ángulo y distancia
    for (int ia = 0; ia <= 0; ++ia) {
        double rot = 135.0 * (1.0 + 0.05*ia);

        for (int id = 0; id <= 0; ++id) {
            double dist = detector_distance_nominal * (1.0 + 0.05*id);

            for (int ix = 0 ; ix <= 0; ++ix) {
                double x_offset = 0.05 * ix; // offset en cm
                for (int iy = 0; iy <=0; ++iy) {
                    double y_offset = 0.05 * iy; // offset en cm
                    GeoConfig cfg;
                    cfg.rot_angle_deg = rot;
                    cfg.det_dist_cm   = dist;
                    cfg.dx_front = x_offset;  // offset en cm
                    cfg.dy_front = y_offset;  // offset en cm
                    cfg.dx_back  = 0.0; // offset en cm
                    cfg.dy_back  = 0.0; // offset en cm

                    cfg.tag = Form("rot%+02d_dist%+02d_dx%+02d_dy%+02d", ia, id, ix, iy);

                    run_simulation(cfg, infile);
                }
            }
        }
    }

    return 0;
}
