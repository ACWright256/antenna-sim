# 2025 0x1ea7bee5
# Copyright @yourmother
# "Ambitiously creating the most obnoxious code since 2002"
# Thanks to:
#   https://www.youtube.com/watch?v=DeG8BaYmvpw
#   https://doi.org/10.1109/8.81458

import numpy as np
import scipy
import pandas
import matplotlib.pyplot as plt

C=299792458

class RotmanSolver:
    """Main container to solve rotman lens parameters :3"""
    def __init__(self, n_beam_ports:int, n_array_ports:int,n_dummy_ports:int, center_f_hz:float):
        """Define basic design constants"""
        self.n_beam_ports = n_beam_ports
        self.n_array_ports = n_array_ports
        self.n_dummy_ports=n_dummy_ports
        self.center_f_hz=center_f_hz
        self.Er=3.55
        

    def update_rotman_params(self,
                      max_steer_ang_deg:float, 
                      focal_angle_deg:float,
                      beta:float,
                      array_elem_spacing:float):
        """Update basic parameters for the rotman lens geometry itself"""

        self.max_steer_ang_rad= np.deg2rad(max_steer_ang_deg)   #theta
        self.focal_angle_rad=np.deg2rad(focal_angle_deg)    #alpha
        self.beta=beta
        
        self.gamma=np.sin(self.max_steer_ang_rad)/np.sin(self.focal_angle_rad) #ratio of sin(theta)/sin(alpha)
        self.lambda0=C/self.center_f_hz                     #target wavelength
        self.f1=2.5*self.lambda0     #focal length based on desired frequency TODO figure out
        self.rho0=1-(1-self.beta**2)/(2*(1-self.beta*np.cos(self.max_steer_ang_rad))) #constant defining the normalized length of a equilateral triangle between two beam ports
        self.array_elem_spacing=array_elem_spacing*self.lambda0 #normalize spacing to wavelength


    def update_trace_paramters(self,
                               trace_width:float,
                               trace_len:float, 
                               beam_taper_len:float,
                               array_taper_len:float,
                               dummy_taper_len:float,
                               polynomial_deg:int,
                               dummy_deviation:float,
                               lens_wid:float
                               ):
        """update trace parameters"""
        self.trace_width=trace_width*np.sqrt(self.Er) #will depend on substrate thickness. normalize to ER
        self.trace_len=trace_len*self.lambda0
        self.beam_taper_len=beam_taper_len*self.lambda0
        self.array_taper_len=array_taper_len*self.lambda0
        self.dummy_taper_len=dummy_taper_len*self.lambda0
        self.polynomial_deg=polynomial_deg #polynomial degree for dummy port approximation 
        self.dummy_deviation=dummy_deviation #deviation of dummy port long axis from lens contour normal
        self.lens_wid=lens_wid #lens width relative to total array of beam ports area width
        self.w0=0.005/self.f1

    def generate_mesh(self):
        """
        This will generate a "mesh" / set of points to be used for importing into kicad :3
        """

    def solve_transmission_lines(self, ws:np.array):
        """
        Given the
        """

    def calculate_normalized_parameters(self):
        """Solve for the basic parameters"""
        self._solve_array_ports()
        self._solve_beam_ports()
  
        plt.figure()
        plt.plot(self.norm_beam_x,self.norm_beam_y,"-o")
        plt.plot(self.norm_array_x,self.norm_array_y, "-o")
        plt.show()
        #print("arrayx ",array_x)
        #print("arrayy ",array_y)
        #print("ws",ws)
        #print("beam_ports_x ",beam_ports_x)
        #print("beam_ports_y ",beam_ports_y)


    def _solve_array_ports(self):
        """solve for the location of the array ports + trace lengths"""
        y0s=self.array_elem_spacing * np.arange(1,self.n_array_ports,1)-(self.n_array_ports)*self.array_elem_spacing/2
        additive = self.n_array_ports%2==0
        if 0.0 not in y0s:
            y0s=np.insert(y0s, self.n_array_ports//2+additive,0)
        zeta=y0s*self.gamma/self.f1
        a=1-((1-self.beta)**2)/((1-self.beta*np.cos(self.focal_angle_rad))**2) - (zeta**2)/self.beta**2
        b=-2*(2*zeta**2)/(self.beta) + 2*(1-self.beta)/(1-self.beta*np.cos(self.focal_angle_rad)) - (zeta**2*(np.sin(self.focal_angle_rad))**2*(1-self.beta))/((1-self.beta*np.cos(self.focal_angle_rad))**2)
        c=-zeta**2 + (zeta**2 * np.sin(self.focal_angle_rad)**2)/(1-self.beta*np.cos(self.focal_angle_rad)) - (zeta**4 * np.sin(self.focal_angle_rad)**4)/(4*(1-self.beta * np.cos(self.focal_angle_rad))**2)
        ws = -(b + np.sqrt(b**2-4*a*c))/(2*a) # length deltas between array port locations and central trace length
        array_y=zeta*(1-ws/self.beta)
        array_x=1-((zeta**2)/2 * np.sin(self.focal_angle_rad)**2 + (1-self.beta)*ws)/(1-self.beta*np.cos(self.focal_angle_rad))
        self.norm_array_x = array_x
        self.norm_array_y = array_y
        self.norm_array_tlen = ws
        #return array_x, array_y, ws #return array port x pos, y pos, and normalized trace lengths

    def _solve_beam_ports(self):
        """solve for location of beam ports"""
        beam_angles = np.linspace(-self.max_steer_ang_rad, self.max_steer_ang_rad, self.n_beam_ports)        
        additive = len(beam_angles)%2==0
        if 0 not in beam_angles:
            beam_angles=np.insert(beam_angles, len(beam_angles)//2+additive,0)
        center_idx = np.where(beam_angles ==0)[0]
        alpha_primes = np.arcsin(np.sin(beam_angles)/self.gamma)    #alpha prime angles for the position of each beam port on a circle
        phis = np.arcsin(np.sin(alpha_primes)*(1-self.rho0)/(self.rho0))      #left over angle
        coordinate_angles = phis+alpha_primes                       
        #initial beam port coordinates from basic geometry
        beam_ports_x=self.rho0*(1-np.cos(coordinate_angles))
        beam_ports_y=self.rho0*np.sin(coordinate_angles)
        
        #beam corrections
        angle_bp_arc_center = np.arctan(beam_ports_y/(self.rho0-beam_ports_x))
        angle_ar_center = np.arctan(beam_ports_y/(1-beam_ports_x))
        
        beam_ports_x_delta = np.array([beam_ports_x[i]-beam_ports_x[i+1] for i in range(len(beam_ports_x)-1)])
        beam_ports_y_delta = np.array([beam_ports_y[i]-beam_ports_y[i+1] for i in range(len(beam_ports_y)-1)])
        angle_bp_arc_delta = np.array([angle_bp_arc_center[i] - angle_bp_arc_center[i+1] for i in range(len(angle_bp_arc_center)-1)])
        
        
        ap_width = (beam_ports_x_delta**2 + beam_ports_y_delta**2) / np.cos(angle_bp_arc_delta)/2
        ap_width=ap_width*self.f1

        Rh = np.sqrt((ap_width+self.trace_width)**2 + self.beam_taper_len**2)
        S = ap_width**2/(8*self.lambda0*Rh)

        corr_factor=self._beam_corr_1(s=S, Rh=Rh)
        offset = len(corr_factor)%2
        corr_factor = np.insert(corr_factor, len(corr_factor)//2+offset,0) #add zero correction to center?
        
        self.norm_beam_x = beam_ports_x + np.cos(angle_ar_center)*corr_factor/self.f1
        self.norm_beam_y = beam_ports_y - np.sin(angle_ar_center)*corr_factor/self.f1

        norm_bp_taper_x =self.norm_beam_x - np.cos(angle_ar_center)*(self.beam_taper_len)/self.f1
        norm_bp_taper_y=self.norm_beam_y + np.sin(angle_ar_center)*(self.beam_taper_len)/self.f1
        norm_bp_taper1_x =self.norm_beam_x - np.cos(angle_ar_center)*(self.trace_len+self.beam_taper_len)/self.f1
        norm_bp_taper1_y=self.norm_beam_y + np.sin(angle_ar_center)*(self.trace_len+self.beam_taper_len)/self.f1

        #self.norm_beam1_x=self.norm_beam_x+np.sin(angle_bp_arc_center)*ap_width/2/self.f1
        #self.norm_beam1_y=self.norm_beam_y+np.sin(angle_bp_arc_center)*ap_width/2/self.f1
        #self.norm_beam2_x=self.norm_beam_x-np.sin(angle_bp_arc_center)*ap_width/2/self.f1
        #self.norm_beam2_y=self.norm_beam_y-np.sin(angle_bp_arc_center)*ap_width/2/self.f1

        self.norm_beam3_x=norm_bp_taper_x+np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam3_y=norm_bp_taper_y+np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam4_x=norm_bp_taper_x-np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam4_y=norm_bp_taper_y-np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1


        self.norm_beam5_x=norm_bp_taper1_x+np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam5_y=norm_bp_taper1_y+np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam6_x=norm_bp_taper1_x-np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1
        self.norm_beam6_y=norm_bp_taper1_y-np.sin(angle_bp_arc_center)*self.trace_width/2/self.f1

    def _beam_corr_1(self, s, Rh):
        """"""
        mask=s>0.06
        print("mask",mask)
        return mask*(2.6746*s**2 + 0.2026*s -0.0085)*Rh
    def _beam_corr_2(self, s, Rh):
        """"""
        mask=s>0.06
        return mask*(-454.9955*s**6 +837.7860*s**5 + 147.3056*s**3 - 9.4030*s**2+0.3223*s)*Rh
