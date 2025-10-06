import numpy as np
import matplotlib.pyplot as plt
from math import degrees

class PPACFissionSimulationSimple:


    def __init__(self, n_events=100000, seed=None, dE_dx_table=None):
        if seed is not None:
            np.random.seed(seed)

        self.n_events = n_events


        self.Z1, self.A1 = 40, 95   # fragment 1: Z=40
        self.Z2, self.A2 = 52, 139  # fragment 2: Z=52
        self.E0_fission = 165.0     # total kinetic energy (MeV)


        self.detector_distance = 5.0   # cm separation between detectors
        self.z_det_back = -self.detector_distance / 2.0
        self.z_det_front = +self.detector_distance / 2.0
        self.detector_width = 20.0     # cm
        self.detector_height = 20.0    # cm
        self.target_radius = 4.        # cm

        self.target_thickness = 1.64e-4   # cm (1.64 µm for uranium target)
        self.backing_thickness = 2.5e-4   # cm (2.5 µm)
        self.gas_thickness = self.detector_distance / 2.0  # cm (2.5 cm)

        if dE_dx_table is None:
            self.dE_dx_table = {
                'Target':       {40: 21.639/1.64*10000 , 52:18.995/1.64*10000 },
                'Backing (Al)': {40: 23.726/2.5*10000 , 52: 19.258/2.5*10000 },
                'Gas (C3F8)':   {40: 3.972/2.5 ,    52: 3.048/2.5 },
            }
        else:
            self.dE_dx_table = dE_dx_table

        self.x_back_local, self.y_back_local = [], []
        self.x_front_local, self.y_front_local = [], []
        self.x_back_world, self.y_back_world = [], []
        self.x_front_world, self.y_front_world = [], []

        self.initial_energy_forward, self.initial_energy_backward = [], []
        self.loss_target_forward, self.loss_target_backward = [], []
        self.loss_backing_forward = []
        self.loss_gas_forward, self.loss_gas_backward = [], []
        self.final_energy_forward, self.final_energy_backward = [], []

        self.theta_deg, self.phi_deg, self.cos_theta = [], [], []


    def target_position(self, n): # uniform in circle, uniform in thickness
        phi = np.random.uniform(0, 2*np.pi, n)
        r = self.target_radius * np.sqrt(np.random.uniform(0, 1, n))
        # random depth inside target
        z_depth = np.random.uniform(-self.target_thickness/2, self.target_thickness/2, n)
        return np.vstack([r*np.cos(phi), r*np.sin(phi), z_depth]).T

    def fission_directions(self, n): # isotropic 
        phi = np.random.uniform(0, 2*np.pi, n)
        theta = np.arccos(np.random.uniform(0, 1, n)) #forward hemisphere
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        return np.vstack([sin_theta*np.cos(phi), sin_theta*np.sin(phi), cos_theta]).T

    def fission_energies(self, n): # based on mass split
        mass_ratio = self.A2 / self.A1
        E1 = self.E0_fission / (1 + mass_ratio)
        E2 = self.E0_fission - E1
        return np.full(n, E1), np.full(n, E2)

    def line_plane_intersection(self, origin, direction, z_plane): #see if the fragment reaches the detector plane
        if abs(direction[2]) < 1e-12:
            return None
        t = (z_plane - origin[2]) / direction[2]
        if t <= 0:
            return None
        return origin + t * direction

    def in_detector(self, point, detector_type='front'): # check if the intersection point is within detector bounds (since they are tilted, x limits depend on the angle by d*tan(45 degrees))
        if detector_type == 'front':
            x_high_limit = self.detector_width / 2 + 2.5
            x_low_limit = -self.detector_width / 2 + 2.5
        else:
            x_high_limit = self.detector_width / 2 - 2.5
            x_low_limit = -self.detector_width / 2 - 2.5
        return (point[0] >= x_low_limit) and (point[0] <= x_high_limit) and (abs(point[1]) <= self.detector_height / 2)

    def primed_to_beam_coordinates(self, point): # rotation by 45 degrees around y axis
        x_p, y_p, z_p = point
        c45 = np.cos(np.radians(45.0))
        s45 = np.sin(np.radians(45.0))
        x_b = x_p * c45 + z_p * s45
        y_b = y_p
        z_b = -x_p * c45 + z_p * s45
        return np.array([x_b, y_b, z_b])

    def energy_loss_simple(self, Z, material_name, path_length_cm):
        if material_name not in self.dE_dx_table:
            return 0.0
        table = self.dE_dx_table[material_name]
        if Z in table:
            dE_dx = table[Z]
        else:
            zs = np.array(list(table.keys()))
            dE_dx = table[int(zs[np.argmin(np.abs(zs - Z))])]
        return dE_dx * path_length_cm

    def simulate(self):
        origins = self.target_position(self.n_events)
        directions = self.fission_directions(self.n_events)
        E1_initial, E2_initial = self.fission_energies(self.n_events)

        accepted = 0
        stopped = 0
        for i in range(self.n_events):
            origin = origins[i]
            direction = directions[i]

            # random assignment of which fragment goes forward
            if np.random.rand() < 0.5:
                forward_frag = {'Z': self.Z1, 'A': self.A1, 'E': E1_initial[i]}
                backward_frag = {'Z': self.Z2, 'A': self.A2, 'E': E2_initial[i]}
            else:
                forward_frag = {'Z': self.Z2, 'A': self.A2, 'E': E2_initial[i]}
                backward_frag = {'Z': self.Z1, 'A': self.A1, 'E': E1_initial[i]}
            dir_forward= direction
            dir_backward= -direction
            cos_f = max(abs(dir_forward[2]), -99999999)

            half_t = self.target_thickness / 2.0
            z0 = origin[2]
            path_target_forward = (half_t - z0) / cos_f
            path_target_backward = (half_t + z0) / cos_f

            # backing: only forward-going fragment
            path_backing = self.backing_thickness / cos_f

            # gas paths
            path_gas_forward = self.gas_thickness / cos_f
            path_gas_backward = self.gas_thickness / cos_f

            loss_target_forward = self.energy_loss_simple(forward_frag['Z'], 'Target', path_target_forward)
            loss_target_backward = self.energy_loss_simple(backward_frag['Z'], 'Target', path_target_backward)

            E_after_target_forward = max(forward_frag['E'] - loss_target_forward, 0.0)
            E_after_target_backward = max(backward_frag['E'] - loss_target_backward, 0.0)
            if E_after_target_forward <= 0.5 or E_after_target_backward <= 0.5:
                stopped += 1
                continue

            loss_back_forward = self.energy_loss_simple(forward_frag['Z'], 'Backing (Al)', path_backing)
            E_after_back_forward = max(E_after_target_forward - loss_back_forward, 0.0)

            E_after_back_backward = E_after_target_backward #backward fragment does not traverse backing

            if E_after_back_forward <= 0.5 or E_after_back_backward <= 0.5:
                stopped += 1
                continue

            # gas energy losses
            loss_gas_forward = self.energy_loss_simple(forward_frag['Z'], 'Gas (C3F8)', path_gas_forward)
            loss_gas_backward = self.energy_loss_simple(backward_frag['Z'], 'Gas (C3F8)', path_gas_backward)

            E_final_forward = max(E_after_back_forward - loss_gas_forward, 0.0)
            E_final_backward = max(E_after_back_backward - loss_gas_backward, 0.0)
            if E_final_forward <= 0.5 or E_final_backward <= 0.5:
                stopped += 1
                continue

            # detector intersections
            hit_front = self.line_plane_intersection(origin, dir_forward, self.z_det_front)
            hit_back = self.line_plane_intersection(origin, dir_backward, self.z_det_back)

            if hit_front is None or hit_back is None:
                continue
            if not (self.in_detector(hit_front, 'front') and self.in_detector(hit_back, 'back')):
                continue

            accepted += 1

            # local positions
            self.x_front_local.append(hit_front[0])
            self.y_front_local.append(hit_front[1])
            self.x_back_local.append(hit_back[0])
            self.y_back_local.append(hit_back[1])

            # beam/world coords
            hf_beam = self.primed_to_beam_coordinates(hit_front)
            hb_beam = self.primed_to_beam_coordinates(hit_back)
            self.x_front_world.append(hf_beam[0])
            self.y_front_world.append(hf_beam[1])
            self.x_back_world.append(hb_beam[0])
            self.y_back_world.append(hb_beam[1])

            

            # store
            self.initial_energy_forward.append(forward_frag['E'])
            self.initial_energy_backward.append(backward_frag['E'])
            self.loss_target_forward.append(loss_target_forward)
            self.loss_target_backward.append(loss_target_backward)
            self.loss_backing_forward.append(loss_back_forward)
            self.loss_gas_forward.append(loss_gas_forward)
            self.loss_gas_backward.append(loss_gas_backward)
            self.final_energy_forward.append(E_final_forward)
            self.final_energy_backward.append(E_final_backward)

            # angles
            frag_vec = hf_beam - hb_beam
            nrm = np.linalg.norm(frag_vec)
            if nrm > 0:
                frag_vec /= nrm
                beam_vec = np.array([0.0, 0.0, 1.0])
                cos_theta = np.dot(frag_vec, beam_vec)
                theta = degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
                phi = np.arctan2(frag_vec[1], frag_vec[0])
                self.cos_theta.append(cos_theta)
                self.theta_deg.append(theta)
                self.phi_deg.append(phi)

        print(f"Simulation completed: {accepted}/{self.n_events} accepted ({accepted/(self.n_events)*100:.2f}%)")
        print(f"Simulation stopped: {stopped}/{self.n_events} stopped ({stopped/self.n_events*100:.2f}%)")
        self._convert_to_arrays()
        return accepted, stopped

    def _convert_to_arrays(self):
        attrs = [
            'x_back_world','y_back_world','x_front_world','y_front_world',
            'x_back_local','y_back_local','x_front_local','y_front_local',
            'theta_deg','phi_deg','cos_theta',
            'initial_energy_forward','initial_energy_backward',
            'loss_target_forward','loss_target_backward','loss_backing_forward',
            'loss_gas_forward','loss_gas_backward','final_energy_forward','final_energy_backward'
        ]
        for a in attrs:
            if isinstance(getattr(self,a), list):
                setattr(self,a,np.array(getattr(self,a)))

    # ---------------- Plotting ----------------
    def plot_detector_hits(self):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        plt.subplots_adjust(wspace=0.3, hspace=0.3)

        def get_limits(x, y, margin=0.5):
            x_min, x_max = np.percentile(x, [0.5, 99.5])
            y_min, y_max = np.percentile(y, [0.5, 99.5])
            return (x_min - margin, x_max + margin), (y_min - margin, y_max + margin)

        # Local coordinates
        xlim, ylim = get_limits(self.x_front_local, self.y_front_local)
        axes[0,0].hist2d(self.x_front_local, self.y_front_local, bins=100, cmap='viridis')
        axes[0,0].set_title('Front detector hits (local)')
        axes[0,0].set_aspect('equal', 'box')
        axes[0,0].set_xlim(xlim)
        axes[0,0].set_ylim(ylim)

        xlim, ylim = get_limits(self.x_back_local, self.y_back_local)
        axes[1,0].hist2d(self.x_back_local, self.y_back_local, bins=100, cmap='viridis')
        axes[1,0].set_title('Back detector hits (local)')
        axes[1,0].set_aspect('equal', 'box')
        axes[1,0].set_xlim(xlim)
        axes[1,0].set_ylim(ylim)

        # World/beam coordinates
        xlim, ylim = get_limits(self.x_front_world, self.y_front_world)
        axes[0,1].hist2d(self.x_front_world, self.y_front_world, bins=100, cmap='viridis')
        axes[0,1].set_title('Front detector hits (beam)')
        axes[0,1].set_aspect('equal', 'box')
        axes[0,1].set_xlim(xlim)
        axes[0,1].set_ylim(ylim)

        xlim, ylim = get_limits(self.x_back_world, self.y_back_world)
        axes[1,1].hist2d(self.x_back_world, self.y_back_world, bins=100, cmap='viridis')
        axes[1,1].set_title('Back detector hits (beam)')
        axes[1,1].set_aspect('equal', 'box')
        axes[1,1].set_xlim(xlim)
        axes[1,1].set_ylim(ylim)



        # 2D histogram of cos(theta) vs phi
        plt.figure(figsize=(7, 5))
        plt.hist2d(self.cos_theta, self.phi_deg, bins=[100, 100], cmap='viridis')
        plt.xlabel('cos(theta)')
        plt.ylabel('phi (rad)')
        plt.colorbar(label='Counts')
        plt.tight_layout()
        plt.show()

    def plot_energy_results(self):
        fig, axes = plt.subplots(2, 1, figsize=(14, 10))

        # 1. Energy loss (forward) vs cos(theta)
        # Distinguish between Z=40 and Z=52 for forward fragment
        Z_forward = np.where(np.array(self.initial_energy_forward) > np.array(self.initial_energy_backward), self.Z1, self.Z2)
        mask_40 = Z_forward == self.Z1
        mask_52 = Z_forward == self.Z2
        # There is no variable x_front_hit defined in this context.
        # If you want to filter by y coordinate of the front detector hit,
        # use self.y_front_world or self.y_front_local as appropriate.
        # For example, to only plot events where y_front_world > 0:

        mask_ypos = abs(np.array(self.x_front_world)) >1e-6

        # Plot for Z=40
        axes[0].scatter(np.array(self.cos_theta)[mask_40 & mask_ypos], np.array(self.loss_target_forward)[mask_40 & mask_ypos], s=0.01, alpha=0.5, label='Target (Z=40)', color='tab:blue')
        axes[0].scatter(np.array(self.cos_theta)[mask_40 & mask_ypos], np.array(self.loss_backing_forward)[mask_40 & mask_ypos], s=0.01, alpha=0.5, label='Backing (Z=40)', color='tab:cyan')
        axes[0].scatter(np.array(self.cos_theta)[mask_40 & mask_ypos], np.array(self.loss_gas_forward)[mask_40 & mask_ypos], s=0.01, alpha=0.5, label='Gas (Z=40)', color='tab:purple')

        # Plot for Z=52
        axes[0].scatter(np.array(self.cos_theta)[mask_52 & mask_ypos], np.array(self.loss_target_forward)[mask_52 & mask_ypos], s=0.01, alpha=0.5, label='Target (Z=52)', color='tab:orange')
        axes[0].scatter(np.array(self.cos_theta)[mask_52 & mask_ypos], np.array(self.loss_backing_forward)[mask_52 & mask_ypos], s=0.01, alpha=0.5, label='Backing (Z=52)', color='tab:olive')
        axes[0].scatter(np.array(self.cos_theta)[mask_52 & mask_ypos], np.array(self.loss_gas_forward)[mask_52 & mask_ypos], s=0.01, alpha=0.5, label='Gas (Z=52)', color='tab:red')

        axes[0].set_xlabel('cos(theta)')
        axes[0].set_ylabel('Energy loss forward (MeV)')
        axes[0].legend(markerscale=10, fontsize=8, ncol=2)

        # 2. Energy loss (backward) vs cos(theta)
        Z_backward = np.where(np.array(self.initial_energy_forward) > np.array(self.initial_energy_backward), self.Z2, self.Z1)
        mask_40_b = Z_backward == self.Z1
        mask_52_b = Z_backward == self.Z2

        # Plot for Z=40
        axes[1].scatter(np.array(self.cos_theta)[mask_40_b & mask_ypos], np.array(self.loss_target_backward)[mask_40_b & mask_ypos], s=0.01, alpha=0.5, label='Target (Z=40)', color='tab:blue')
        axes[1].scatter(np.array(self.cos_theta)[mask_40_b & mask_ypos], np.array(self.loss_gas_backward)[mask_40_b & mask_ypos], s=0.01, alpha=0.5, label='Gas (Z=40)', color='tab:purple')

        # Plot for Z=52
        axes[1].scatter(np.array(self.cos_theta)[mask_52_b & mask_ypos], np.array(self.loss_target_backward)[mask_52_b & mask_ypos], s=0.01, alpha=0.5, label='Target (Z=52)', color='tab:orange')
        axes[1].scatter(np.array(self.cos_theta)[mask_52_b & mask_ypos], np.array(self.loss_gas_backward)[mask_52_b & mask_ypos], s=0.01, alpha=0.5, label='Gas (Z=52)', color='tab:red')

        axes[1].set_xlabel('cos(theta)')
        axes[1].set_ylabel('Energy loss backward (MeV)')
        axes[1].legend(markerscale=10, fontsize=8, ncol=2)
        plt.tight_layout()
        plt.show()

        # 5. Plot total energy loss and its contributions as a function of Z (charge) for forward and backward fragments
        # as histograms overlapped using step

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # 3 & 4. Histograms of energy losses by Z
        colors = plt.get_cmap('tab10').colors
        bins = 200

        # Reconstruct which Z was forward/backward for each event
        Z_forward = np.where(np.array(self.initial_energy_forward) > np.array(self.initial_energy_backward), self.Z1, self.Z2)
        Z_backward = np.where(Z_forward == self.Z1, self.Z2, self.Z1)

        # Forward fragment: total energy loss (target + gas + backing), by Z
        loss_forward = np.array(self.loss_target_forward) + np.array(self.loss_gas_forward) + np.array(self.loss_backing_forward)
        # Backward fragment: total energy loss (target + gas), by Z
        loss_backward = np.array(self.loss_target_backward) + np.array(self.loss_gas_backward)

        # Find global min/max for bin range
        min_loss = min(loss_forward.min(), loss_backward.min())
        max_loss = max(loss_forward.max(), loss_backward.max())
        bin_range = (min_loss, max_loss)

        # Forward fragment
        mask_40 = Z_forward == self.Z1
        mask_52 = Z_forward == self.Z2

        # For each Z, plot histograms of each contribution and total
        for z, mask, color in zip([self.Z1, self.Z2], [mask_40, mask_52], ['tab:blue', 'tab:orange']):
            axes[0].hist(np.array(self.loss_target_forward)[mask], bins=100, histtype='step', label=f'Target (Z={z})', color=color, linewidth=1.2, alpha=0.8)
            axes[0].hist(np.array(self.loss_backing_forward)[mask], bins=100, histtype='step', label=f'Backing (Z={z})', color=color, linestyle='--', linewidth=1.2, alpha=0.8)
            axes[0].hist(np.array(self.loss_gas_forward)[mask], bins=100, histtype='step', label=f'Gas (Z={z})', color=color, linestyle=':', linewidth=1.2, alpha=0.8)
            axes[0].hist(loss_forward[mask], bins=100, histtype='step', label=f'Total (Z={z})', color=color, linestyle='-.', linewidth=1.5, alpha=1.0)

        axes[0].set_ylabel('Counts')
        axes[0].legend(fontsize=8, ncol=2)

        # Backward fragment
        
        mask_40_b = Z_backward == self.Z1
        mask_52_b = Z_backward == self.Z2

        for z, mask, color in zip([self.Z1, self.Z2], [mask_40_b, mask_52_b], ['tab:blue', 'tab:orange']):
            axes[1].hist(np.array(self.loss_target_backward)[mask], bins=100, histtype='step', label=f'Target (Z={z})', color=color, linewidth=1.2, alpha=0.8)
            axes[1].hist(np.array(self.loss_gas_backward)[mask], bins=100, histtype='step', label=f'Gas (Z={z})', color=color, linestyle=':', linewidth=1.2, alpha=0.8)
            axes[1].hist(loss_backward[mask], bins=100, histtype='step', label=f'Total (Z={z})', color=color, linestyle='-.', linewidth=1.5, alpha=1.0)

        axes[1].set_xlabel('Energy loss (MeV)')
        axes[1].set_ylabel('Counts')
        axes[1].legend(fontsize=8, ncol=2)

        plt.tight_layout()
        plt.show()

    def print_summary(self):
        print("\n"+"="*60)
        print("SIMULATION SUMMARY (Z=40 and Z=52 fragments)")
        print("="*60)
        n_acc = len(self.final_energy_forward)
        print(f"Accepted events: {n_acc}/{self.n_events} ({n_acc/self.n_events*100:.2f}%)")


# ---------------- Example usage ----------------
if __name__ == "__main__":
    sim = PPACFissionSimulationSimple(n_events=2000000, seed=42)
    sim.simulate()
    sim.plot_detector_hits()
    sim.plot_energy_results()
    sim.print_summary()
