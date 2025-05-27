from rotman_solver import RotmanSolver

def main():
    """"""
    rotman_solver=RotmanSolver(5,7,3,10.25e9)
    #rotman_solver.update_rotman_params(max_steer_ang_deg=30, focal_angle_deg=30, beta=0.9, array_elem_spacing=0.5)
    rotman_solver.update_rotman_params(max_steer_ang_deg=35, focal_angle_deg=30, beta=0.9, array_elem_spacing=0.5)
    rotman_solver.update_trace_paramters(trace_width=0.7e-3, trace_len=1,beam_taper_len=3, array_taper_len=3, dummy_taper_len=3, polynomial_deg=2, dummy_deviation=0.8,lens_wid=1.4)
    rotman_solver.calculate_normalized_parameters()
if __name__ == '__main__':
    main()