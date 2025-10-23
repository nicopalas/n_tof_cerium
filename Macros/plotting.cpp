#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <algorithm>
#include <vector>

void make_histograms() {
    gStyle->SetOptStat(0);

    // Input / output
    TFile *fin = new TFile("fission_events_cerium_symmetric.root", "READ");
    TTree *t = (TTree*)fin->Get("fission_process");
    TFile *fout = new TFile("fission_histograms_cerium_symmetric.root", "RECREATE");

    // Helper lambda to get range of variable
    auto get_range = [&](const char* var){
        t->Draw(Form("%s>>htemp(1000)", var), "", "goff");
        TH1 *h = (TH1*)gDirectory->Get("htemp");
        double min = h->GetXaxis()->GetXmin();
        double max = h->GetXaxis()->GetXmax();
        delete h;
        return std::pair<double,double>(min, max);
    };

    // =====================================================================
    // (1) Forward fragment losses (target, backing, gas)
    // =====================================================================
    auto r_target_f = get_range("loss_target_forward");
    auto r_backing_f = get_range("loss_backing_forward");
    auto r_gas_f = get_range("loss_gas_forward");

    double xmin_f = std::min({r_target_f.first, r_backing_f.first, r_gas_f.first});
    double xmax_f = std::max({r_target_f.second, r_backing_f.second, r_gas_f.second});

    TH1F *h_target_f = new TH1F("h_target_f",
        "#font[12]{Forward Fragment Energy Losses};#DeltaE (MeV);Counts",150,xmin_f,xmax_f);
    TH1F *h_backing_f = new TH1F("h_backing_f","",150,xmin_f,xmax_f);
    TH1F *h_gas_f     = new TH1F("h_gas_f","",150,xmin_f,xmax_f);

    t->Draw("loss_target_forward>>h_target_f","","goff");
    t->Draw("loss_backing_forward>>h_backing_f","","goff");
    t->Draw("loss_gas_forward>>h_gas_f","","goff");

    // compute true global max
    double max_target_f = h_target_f->GetBinContent(h_target_f->GetMaximumBin());
    double max_backing_f = h_backing_f->GetBinContent(h_backing_f->GetMaximumBin());
    double max_gas_f     = h_gas_f->GetBinContent(h_gas_f->GetMaximumBin());
    double ymax_f = std::max({max_target_f, max_backing_f, max_gas_f}) * 1.1;

    // Style
    h_target_f->SetLineColor(kRed);
    h_backing_f->SetLineColor(kBlue);
    h_gas_f->SetLineColor(kGreen+2);
    h_target_f->SetLineWidth(2);
    h_backing_f->SetLineWidth(2);
    h_gas_f->SetLineWidth(2);
    h_target_f->SetMaximum(ymax_f);
    h_backing_f->SetMaximum(ymax_f);
    h_gas_f->SetMaximum(ymax_f);

    // Draw
    TCanvas *c1 = new TCanvas("c1","Forward Losses",800,600);
    h_target_f->Draw("HIST");
    h_backing_f->Draw("HIST SAME");
    h_gas_f->Draw("HIST SAME");

    auto leg1 = new TLegend(0.65,0.7,0.88,0.88);
    leg1->SetHeader("#font[12]{Forward losses}","C");
    leg1->AddEntry(h_target_f,"Target","l");
    leg1->AddEntry(h_backing_f,"Backing","l");
    leg1->AddEntry(h_gas_f,"Gas","l");
    leg1->Draw();
    c1->Write();

    // =====================================================================
    // (2) Backward fragment losses (target, gas)
    // =====================================================================
    auto r_target_b = get_range("loss_target_backward");
    auto r_gas_b = get_range("loss_gas_backward");

    double xmin_b = std::min(r_target_b.first, r_gas_b.first);
    double xmax_b = std::max(r_target_b.second, r_gas_b.second);

    TH1F *h_target_b = new TH1F("h_target_b",
        "#font[12]{Backward Fragment Energy Losses};#DeltaE (MeV);Counts",150,xmin_b,xmax_b);
    TH1F *h_gas_b = new TH1F("h_gas_b","",150,xmin_b,xmax_b);

    t->Draw("loss_target_backward>>h_target_b","","goff");
    t->Draw("loss_gas_backward>>h_gas_b","","goff");

    double max_target_b = h_target_b->GetBinContent(h_target_b->GetMaximumBin());
    double max_gas_b    = h_gas_b->GetBinContent(h_gas_b->GetMaximumBin());
    double ymax_b = std::max(max_target_b, max_gas_b) * 1.1;

    h_target_b->SetLineColor(kRed);
    h_gas_b->SetLineColor(kGreen+2);
    h_target_b->SetLineWidth(2);
    h_gas_b->SetLineWidth(2);
    h_target_b->SetMaximum(ymax_b);
    h_gas_b->SetMaximum(ymax_b);

    TCanvas *c2 = new TCanvas("c2","Backward Losses",800,600);
    h_target_b->Draw("HIST");
    h_gas_b->Draw("HIST SAME");

    auto leg2 = new TLegend(0.65,0.75,0.88,0.88);
    leg2->SetHeader("#font[12]{Backward losses}","C");
    leg2->AddEntry(h_target_b,"Target","l");
    leg2->AddEntry(h_gas_b,"Gas","l");
    leg2->Draw();
    c2->Write();

    // =====================================================================
    // (3) Gas losses for cathode_x+cathode_y (forward/backward)
    // =====================================================================
    t->Draw("(loss_gas_cathode_x_f+loss_gas_cathode_y_f)>>h_temp_f(150)","","goff");
    t->Draw("(loss_gas_cathode_x_b+loss_gas_cathode_y_b)>>h_temp_b(150)","","goff");
    TH1F *hf = (TH1F*)gDirectory->Get("h_temp_f");
    TH1F *hb = (TH1F*)gDirectory->Get("h_temp_b");

    double xmin_c = std::min(hf->GetXaxis()->GetXmin(), hb->GetXaxis()->GetXmin());
    double xmax_c = std::max(hf->GetXaxis()->GetXmax(), hb->GetXaxis()->GetXmax());

    TH1F *h_gas_f_cath = new TH1F("h_gas_f_cath",
        "#font[12]{Gas Loss Cathodes};#DeltaE_{cath} (MeV);Counts",150,xmin_c,xmax_c);
    TH1F *h_gas_b_cath = new TH1F("h_gas_b_cath","",150,xmin_c,xmax_c);
    t->Draw("(loss_gas_cathode_x_f+loss_gas_cathode_y_f)>>h_gas_f_cath","","goff");
    t->Draw("(loss_gas_cathode_x_b+loss_gas_cathode_y_b)>>h_gas_b_cath","","goff");

    double max_gas_f_cath = h_gas_f_cath->GetBinContent(h_gas_f_cath->GetMaximumBin());
    double max_gas_b_cath = h_gas_b_cath->GetBinContent(h_gas_b_cath->GetMaximumBin());
    double ymax_c = std::max(max_gas_f_cath, max_gas_b_cath) * 1.3;

    h_gas_f_cath->SetLineColor(kBlue);
    h_gas_b_cath->SetLineColor(kRed);
    h_gas_f_cath->SetLineWidth(2);
    h_gas_b_cath->SetLineWidth(2);
    h_gas_f_cath->SetMaximum(ymax_c);
    h_gas_b_cath->SetMaximum(ymax_c);

    TCanvas *c3 = new TCanvas("c3","Gas Loss Cathodes",800,600);
    h_gas_f_cath->Draw("HIST");
    h_gas_b_cath->Draw("HIST SAME");
    auto leg3 = new TLegend(0.65,0.75,0.88,0.88);
    leg3->AddEntry(h_gas_f_cath,"Forward","l");
    leg3->AddEntry(h_gas_b_cath,"Backward","l");
    leg3->Draw();
    c3->Write();

    // =====================================================================
    // Remaining plots (4–13): single or 2D, unchanged
    // =====================================================================
    TCanvas *c4 = new TCanvas("c4","TOF vs Cathode",800,600);
    t->Draw("tof_diff:(loss_gas_cathode_x_f+loss_gas_cathode_y_f)>>h_tof_vs_cath(150,0,0,150,-20,20)",
            "tof_diff>-20 && tof_diff<20","COLZ");
    ((TH2F*)gDirectory->Get("h_tof_vs_cath"))
        ->SetTitle("#font[12]{t_{1}-t_{0} vs Cathode Sum (Forward)};#DeltaE_{cath} (MeV);t_{1}-t_{0} (ns)");
    c4->Write();

    TCanvas *c5 = new TCanvas("c5","tof_diff",800,600);
    t->Draw("tof_diff>>h_tof(150,-20,20)","tof_diff>-20 && tof_diff<20");
    ((TH1F*)gDirectory->Get("h_tof"))
        ->SetTitle("#font[12]{Time of Flight Difference};t_{1}-t_{0} (ns);Counts");
    c5->Write();

    TCanvas *c6 = new TCanvas("c6","tof_diff vs theta",800,600);
    t->Draw("tof_diff:theta_deg>>h_tof_vs_theta(150,0,0,150,-20,20)",
            "tof_diff>-20 && tof_diff<20","COLZ");
    ((TH2F*)gDirectory->Get("h_tof_vs_theta"))
        ->SetTitle("#font[12]{t_{1}-t_{0} vs #theta};#theta (^{#circ});t_{1}-t_{0} (ns)");
    c6->Write();

    TCanvas *c7 = new TCanvas("c7","phi vs theta",800,600);
    t->Draw("phi_deg:theta_deg>>h_phi_vs_theta(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_phi_vs_theta"))
        ->SetTitle("#font[12]{#phi vs #theta};#theta (^{#circ});#phi (^{#circ})");
    c7->Write();

    TCanvas *c8 = new TCanvas("c8","Local coordinates",1000,400);
    c8->Divide(2,1);
    c8->cd(1);
    t->Draw("y_back_local:x_back_local>>h_back_local(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_back_local"))
        ->SetTitle("#font[12]{Backward Local};x_{back} (cm);y_{back} (cm)");
    c8->cd(2);
    t->Draw("y_front_local:x_front_local>>h_front_local(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_front_local"))
        ->SetTitle("#font[12]{Forward Local};x_{front} (cm);y_{front} (cm)");
    c8->Write();

    TCanvas *c9 = new TCanvas("c9","World coordinates",1000,400);
    c9->Divide(2,1);
    c9->cd(1);
    t->Draw("y_back_world:x_back_world>>h_back_world(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_back_world"))
        ->SetTitle("#font[12]{Backward World};x_{back} (cm);y_{back} (cm)");
    c9->cd(2);
    t->Draw("y_front_world:x_front_world>>h_front_world(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_front_world"))
        ->SetTitle("#font[12]{Forward World};x_{front} (cm);y_{front} (cm)");
    c9->Write();

    TCanvas *c10 = new TCanvas("c10","Energy Backward",800,600);
    t->Draw("final_energy_backward_mevu:initial_energy_backward_mevu>>h_energy_b(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_energy_b"))
        ->SetTitle("#font[12]{Initial vs Final Energy (Backward)};E_{i} (MeV/u);E_{f} (MeV/u)");
    c10->Write();

    TCanvas *c11 = new TCanvas("c11","Energy Forward",800,600);
    t->Draw("final_energy_forward_mevu:initial_energy_forward_mevu>>h_energy_f(150,0,0,150,0,0)","","COLZ");
    ((TH2F*)gDirectory->Get("h_energy_f"))
        ->SetTitle("#font[12]{Initial vs Final Energy (Forward)};E_{i} (MeV/u);E_{f} (MeV/u)");
    c11->Write();

    // =====================================================================
    // NEW: Time of flight difference vs forward fragment energies
    // =====================================================================
    
    // TOF difference vs initial energy of forward fragment
    TCanvas *c12 = new TCanvas("c12","TOF vs Initial Energy Forward",800,600);
    t->Draw("tof_diff:initial_energy_forward_mevu>>h_tof_vs_ini_energy_f(150,0,0,150,-20,20)",
            "tof_diff>-20 && tof_diff<20","COLZ");
    ((TH2F*)gDirectory->Get("h_tof_vs_ini_energy_f"))
        ->SetTitle("#font[12]{t_{1}-t_{0} vs Initial Energy (Forward)};E_{i} (MeV/u);t_{1}-t_{0} (ns)");
    c12->Write();

    // TOF difference vs final energy of forward fragment  
    TCanvas *c13 = new TCanvas("c13","TOF vs Final Energy Forward",800,600);
    t->Draw("tof_diff:final_energy_forward_mevu>>h_tof_vs_fin_energy_f(150,0,0,150,-20,20)",
            "tof_diff>-20 && tof_diff<20","COLZ");
    ((TH2F*)gDirectory->Get("h_tof_vs_fin_energy_f"))
        ->SetTitle("#font[12]{t_{1}-t_{0} vs Final Energy (Forward)};E_{f} (MeV/u);t_{1}-t_{0} (ns)");
    c13->Write();

    // ---------------------------------------------------------------------
        TCanvas *c14 = new TCanvas("c14", "in_Z1 and in_Z2", 800, 600);
        TH1F *h_in_Z1 = new TH1F("h_in_Z1", "#font[12]{Z_{back} and Z_{forward}};Z;Counts", 100, 0, 0);
        TH1F *h_in_Z2 = new TH1F("h_in_Z2", "", 100, 0, 0);

        t->Draw("backward_Z>>h_in_Z1", "abs(tof_diff)<1e9", "goff");
        t->Draw("forward_Z>>h_in_Z2", "abs(tof_diff)<1e9", "goff");

        h_in_Z1->SetLineColor(kBlue+1);
        h_in_Z2->SetLineColor(kRed+1);
        h_in_Z1->SetLineWidth(2);
        h_in_Z2->SetLineWidth(2);

        double max_in_Z1 = h_in_Z1->GetBinContent(h_in_Z1->GetMaximumBin());
        double max_in_Z2 = h_in_Z2->GetBinContent(h_in_Z2->GetMaximumBin());
        double ymax_z = std::max(max_in_Z1, max_in_Z2) * 1.2;
        h_in_Z1->SetMaximum(ymax_z);
        h_in_Z2->SetMaximum(ymax_z);

        h_in_Z1->Draw("HIST");
        h_in_Z2->Draw("HIST SAME");
        h_in_Z1->GetXaxis()->SetLimits(15, 48);
        h_in_Z2->GetXaxis()->SetLimits(15, 48);

        auto leg_z = new TLegend(0.7, 0.75, 0.88, 0.88);
        leg_z->AddEntry(h_in_Z1, "Z_{back}", "l");
        leg_z->AddEntry(h_in_Z2, "Z_{forward}", "l");
        leg_z->Draw();

        c14->Write();
    fout->Write();
    fout->Close();
    fin->Close();
}
