
import os
import numpy as np
from openmdao.main.api import VariableTree, Container, Component
from openmdao.lib.datatypes.api import Int, Str, Float, List, Array, Enum, Bool, VarTree, Dict
from vartrees import *
from hawc2_vartrees import *

import sys


def _get_fmt(val):

    def _get_fmt1(val):
        if isinstance(val, str):
            return '%s'
        elif isinstance(val, int):
            return '%i'
        elif isinstance(val, float):
            return '%.16f'

    if isinstance(val, list):
        return ' '.join([_get_fmt1(v) for v in val])
    else:
        return _get_fmt1(val)

def write_pcfile(path, pc):

    fid = open(path, 'w')
    fid.write('%i %s\n' % (pc.nset, pc.desc))
    for i in range(pc.nset):
        pcset = pc.pc_sets[i]
        fid.write('%i\n' % pcset.np)
        for j in range(pcset.np):
            polar = pcset.polars[j]
            fid.write('%i %i %f %s\n' % (j + 1, polar.aoa.shape[0], polar.rthick, polar.desc))
            for k in range(polar.aoa.shape[0]):
                fid.write('% 20.12e  % 20.12e  % 20.12e  % 20.12e \n' % \
                          (polar.aoa[k], polar.cl[k], polar.cd[k], polar.cm[k]))
    fid.close()

def write_aefile(path, b):
    """write the blade shape to the ae_filename"""

    fid = open(path, 'w')

    fid.write('1 HAWC_AE data\n')
    fid.write("1 %i r,c,t/c prin.set\n" % (b.s.shape[0]))
    data = np.array([b.s,
                     b.chord,
                     np.minimum(100., b.rthick),
                     b.aeset]).T
    np.savetxt(fid, data, fmt="%f %f %f %i")
    fid.close()

def write_stfile(path, body, case_id):

    """write the beam structural data to an st_filename"""
    if body.st_input_type is 0:
        header = ['r', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh', 'E',
                  'G', 'I_x', 'I_y', 'K', 'k_x', 'k_y', 'A', 'pitch', 'x_e', 'y_e']
        # for readable files with headers above the actual data column
        exp_prec = 10             # exponential precesion
        col_width = exp_prec + 8  # column width required for exp precision
        header_full = '='*20*col_width + '\n'
        header_full += ''.join([(hh + ' [%i]').center(col_width+1)%i for i, hh in enumerate(header)])+'\n'
        header_full += '='*20*col_width + '\n'
    else:
        header = ['r', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_e', 'y_e', 'K_11', 
                  'K_12', 'K_13', 'K_14', 'K_15', 'K_16', 'K_22', 'K_23', 
                  'K_24', 'K_25', 'K_26', 'K_33', 'K_34', 'K_35', 'K_36', 
                  'K_44', 'K_45', 'K_46',
                  'K_55', 'K_56', 'K_66']
        # for readable files with headers above the actual data column
        exp_prec = 10             # exponential precesion
        col_width = exp_prec + 8  # column width required for exp precision
        header_full = '='*31*col_width + '\n'
        header_full += ''.join([(hh + ' [%i]').center(col_width+1)%i for i, hh in enumerate(header)])+'\n'
        header_full += '='*31*col_width + '\n'

    fid = open(path, 'w')
    fid.write('1  number of sets, Nset\n' % body.body_set[1])
    fid.write('-----------------\n')
    fid.write('#1 written using the HAWC2 OpenMDAO wrapper\n')
    fid.write('Case ID: %s\n' % case_id)
    if body.st_input_type is 0:
    	for i, st in enumerate(body.beam_structure):
        	fid.write(header_full)
        	fid.write('$%i %i\n' % (i + 1, st.s.shape[0]))
        	data = np.array([st.s,
                          st.dm,
                	         st.x_cg,
                          st.y_cg,
                          st.ri_x,
                          st.ri_y,
    	                    st.x_sh,
                          st.y_sh,
                          st.E,
                	         st.G,
                          st.I_x,
                          st.I_y,
                          st.K,
                          st.k_x,
                          st.k_y,
                          st.A,
                	         st.pitch,
                          st.x_e,
                          st.y_e]).T
        	np.savetxt(fid, data, fmt='%'+' %i.%ie' % (col_width, exp_prec) )
    else:
        for i, st in enumerate(body.beam_structure):
            fid.write(header_full)
            fid.write('$%i %i\n' % (i + 1, st.s.shape[0]))
            data = np.array([st.s,
                             st.dm,
                             st.x_cg,
                             st.y_cg,
                             st.ri_x,
                             st.ri_y,
                             st.x_e,
                             st.y_e,
                             st.K_11,
                             st.K_12,
                             st.K_13,
                             st.K_14,
                             st.K_15,
                             st.K_16,
                             st.K_22,
                             st.K_23,
                             st.K_24,
                             st.K_25,
                             st.K_26,
                             st.K_33,
                             st.K_34,
                             st.K_35,
                             st.K_36,
                             st.K_44,
                             st.K_45,
                             st.K_46,
                             st.K_55,
                             st.K_56,
                             st.K_66]).T
            np.savetxt(fid, data, fmt='%'+' %i.%ie' % (col_width, exp_prec) )
	fid.close()

class HAWC2InputWriter(Component):

    from_file = Bool(True, iotype='in')
    case_id = Str('hawc2_case', iotype='in')
    htc_master = List(iotype='out', desc='Final master list of inputs')

    data_directory = Str('data', iotype='in')
    res_directory = Str(iotype='in')
    turb_directory = Str(iotype='in')
    log_directory = Str(iotype='in')
    control_directory = Str('control', iotype='in')

    vartrees = VarTree(HAWC2VarTrees(), iotype='in')
    Flag = Bool(False, iotype='out')


    def __init__(self):
        super(HAWC2InputWriter, self).__init__()

        self.force_execute = True

        self.master = []
        self.structure = []
        self.aerodrag = []
        self.controlinp = []
        self.sensors = []

        self._nbodies = 0
        self.basedir = os.getcwd()                        ### tkba: addition for flaps ###
        if not os.path.exists(self.data_directory):
            os.mkdir(self.data_directory)


    def write_master(self):

        #print 'writing case id %s' % self.case_id
        fid = open(self.case_id + '.htc', 'w')

        for i, line in enumerate(self.master):
                line = line.rstrip()+'\n'
                self.master[i] = line
                fid.write(line)
        fid.close()

        # too lazy to change it
        self.htc_master = self.master


    def write_all(self):

        self.master = []
        self.structure = []
        self.aerodrag = []
        self.controlinp = []
        self.sensors = []

        self.configure_wt()
        self.write_simulation()
        self.write_wind()
        self.write_aero()
        self.write_aerodrag()
        self.write_structure_res()
        self.update_c12axis()
        self.write_main_bodies()
        self.write_orientations()
        self.write_constraints()
        self.structure.insert(0, 'begin new_htc_structure;')
        self.structure.append('end new_htc_structure;')
        self.write_control_dll()
        if self.from_file:
            self.write_output_fromfile()
        else:
            self.write_output()
        #self.write_master()

        self.master.extend(self.aerodrag)
        self.master.extend(self.structure)
        self.master.extend(self.controlinp)
        self.master.extend(self.sensors)

    def execute(self):

        self._logger.info('writing HAWC2 input files to disk ...')

        if not os.path.exists(self.data_directory):
            os.mkdir(self.data_directory)
        self.case_idout = self.case_id
        self.write_all()
        self.write_master()
        self.write_pcfile()
        self.write_aefile()
        self.Flag = True

    def configure_wt(self):
        if not self.from_file:
            self.configure_tower_body()
            self.configure_towertop_body()
            self.configure_shaft_body()
            self.configure_hub_bodies()
            self.configure_blade_bodies()

            self.add_tower_aerodrag()
            self.add_nacelle_aerodrag()

    def add_tower_aerodrag(self, cd=0.6):
        """convenience function to add tower drag"""
        geom = np.zeros((2, 2))
        geom[0, :] = 0, self.vartrees.tower.bottom_diameter
        geom[1, :] = self.vartrees.tower.height, self.vartrees.tower.top_diameter
        self.add_aerodrag_element('tower', geom, cd)

    def add_nacelle_aerodrag(self, cd=0.8, width=10.):
        """convenience function to add nacelle drag"""

        geom = np.zeros((2, 2))

        shaft = self.vartrees.main_bodies.get_main_body('shaft')

        geom[0, :] = 0, shaft.c12axis[0, 2]
        geom[1, :] = shaft.c12axis[-1, 2], width
        self.add_aerodrag_element('shaft', geom, cd)

    def configure_tower_body(self):
        """convenience method for adding tower body with orientation and constraints"""
        b = self.vartrees.main_bodies.add_main_body('tower')
        b.c12axis = np.zeros((10, 4))
        b.c12axis[:, 2] = np.linspace(0, -self.tower.height, 10)
        b.add_orientation('base')
        b.orientations[0].eulerang.append(np.array([0, 0, 0]))
        b.add_constraint('fixed')
        print 'Not sure it makes sense configure_tower_body'
        return b

    def configure_towertop_body(self):
        """convenience method for adding towertop body with orientation and constraints"""

        b = self.vartrees.main_bodies.add_main_body('towertop')
        b.c12axis = np.zeros((2, 4))
        b.c12axis[-1, 2] = -self.vartrees.nacelle.diameter / 2.
        b.add_orientation('relative')
        b.orientations[0].mbdy1_name = 'tower'
        b.orientations[0].eulerang.append(np.array([0, 0, 0]))
        b.add_constraint('fixed_to_body', body1='tower')
        print 'Not sure it makes sense configure_towertop_body'

    def configure_shaft_body(self):
        """convenience method for adding shaft body with orientation and constraints"""

        b =self.vartrees.main_bodies.add_main_body('shaft')
        b.c12axis = np.zeros((5, 4))
        b.c12axis[:, 2] = np.linspace(0, self.vartrees.shaft.length, 5)
        b.add_orientation('relative')
        b.orientations[0].mbdy1_name = 'towertop'
        b.orientations[0].eulerang.append(np.array([90, 0, 0]))
        b.orientations[0].eulerang.append(np.array([self.rotor.tilt_angle, 0, 0]))
        b.orientations[0].mbdy2_ini_rotvec_d1 = 0.314
        b.orientations[0].rotation_dof = [0, 0, -1]
        b.add_constraint('free', body1='towertop', con_name='shaft_rot', DOF=np.array([0,0,0,0,0,-1]))

    def configure_hub_bodies(self):
        """convenience method for adding hub bodies with orientation and constraints"""

        b = self.vartrees.main_bodies.add_main_body('hub1')
        b.c12axis = np.zeros((2, 4))
        b.c12axis[1, 2] = self.vartrees.hub.diameter/2.
        b.nbodies = 1
        b.add_orientation('relative')
        b.orientations[0].mbdy1_name = 'shaft'
        b.orientations[0].eulerang.append(np.array([-90, 0, 0]))
        b.orientations[0].eulerang.append(np.array([0., 180., 0]))
        b.orientations[0].eulerang.append(np.array([self.vartrees.rotor.cone_angle, 0, 0]))
        b.add_constraint('fixed_to_body', body1='shaft')

        for i in range(1, self.vartrees.rotor.nblades):
            b = self.vartrees.main_bodies.add_main_body('hub'+str(i+1))
            b.copy_main_body = 'hub1'
            b.add_orientation('relative')
            b.orientations[0].mbdy1_name = 'shaft'
            b.orientations[0].eulerang.append(np.array([-90, 0, 0]))
            b.orientations[0].eulerang.append(np.array([0., 60. - (i-1) * 120., 0]))
            b.orientations[0].eulerang.append(np.array([self.vartrees.rotor.cone_angle, 0, 0]))
            b.add_constraint('fixed_to_body', body1='shaft')

    def configure_blade_bodies(self):
        """convenience method for adding blade bodies with orientation and constraints"""

        b = self.vartrees.main_bodies.add_main_body('blade1')
        b.c12axis[:, :3] = self.vartrees.blade_ae.c12axis
        b.c12axis[:,  3] = self.vartrees.blade_ae.twist
        b.nbodies = 10
        b.add_orientation('relative')
        b.orientations[0].mbdy1_name = 'hub1'
        b.orientations[0].eulerang.append(np.array([0, 0, 0]))
        b.add_constraint('prescribed_angle', body1='hub1', con_name='pitch1', DOF=np.array([0,0,0,0,0,-1]))

        for i in range(1, self.rotor.nblades):
            b = self.vartrees.main_bodies.add_main_body('blade'+str(i+1))
            b.copy_main_body = 'blade1'
            b.add_orientation('relative')
            b.orientations[0].mbdy1_name = 'hub'+str(i+1)
            b.orientations[0].eulerang.append(np.array([0, 0, 0]))
            b.add_constraint('prescribed_angle', body1='hub'+str(i+1), con_name='pitch'+str(i+1), DOF=np.array([0,0,0,0,0,-1]))

    def write_simulation(self):

        sim = []
        sim.append('begin simulation ;')
        sim.append('  time_stop    %s ;' % self.vartrees.sim.time_stop)
        sim.append('  solvertype   %i ;    (newmark)' % self.vartrees.sim.solvertype)
        # sim.append(';  animation %s ;' % (os.path.join(self.res_directory, self.case_id+'_animation.dat')))
        sim.append('  convergence_limits %3.6f %3.6f %3.6f;' % (self.vartrees.sim.convergence_limits[0],
                                                                self.vartrees.sim.convergence_limits[1],
                                                                self.vartrees.sim.convergence_limits[2]))
        sim.append('  on_no_convergence continue ;')
        sim.append('  max_iterations %i ;' % self.vartrees.sim.max_iterations)
        sim.append('  logfile %s ;' % (os.path.join(self.log_directory,self.case_id+'.log')))
        sim.append('  begin newmark ;')
        sim.append('    deltat    %1.3f ;' % self.vartrees.sim.newmark_deltat)
        sim.append('  end newmark ;')
        sim.append('end simulation ;')

        self.master.extend(sim)

    def write_wind(self):

        wind = []
        wind.append('begin wind ;')
        wind.append('  density                %3.6f ;' % self.vartrees.wind.density)
        wind.append('  wsp                    %3.6f ;' % self.vartrees.wind.wsp)
        wind.append('  tint                   %3.6f ;' % self.vartrees.wind.tint)
        wind.append('  horizontal_input       %i ; 0=false, 1=true' %
                    self.vartrees.wind.horizontal_input)
        wind.append('  windfield_rotations    %3.6f %3.6f %3.6f ;' %
                    (self.vartrees.wind.windfield_rotations[0],
                     self.vartrees.wind.windfield_rotations[1],
                     self.vartrees.wind.windfield_rotations[2]))
        wind.append('  center_pos0            %3.6f %3.6f %3.6f ; hub height' %
                    (self.vartrees.wind.center_pos0[0],
                     self.vartrees.wind.center_pos0[1],
                     self.vartrees.wind.center_pos0[2]))
        wind.append('  shear_format           %i %3.6f ;' %
                    (self.vartrees.wind.shear_type,
                     self.vartrees.wind.shear_factor))
        wind.append('  turb_format            %i ;  0=none, 1=mann, 2=flex' %
                    self.vartrees.wind.turb_format)
        wind.append('  tower_shadow_method    %i ;  0=none, 1=potential flow, \
                    2=jet' % self.vartrees.wind.tower_shadow_method)

        for wind_ramp in self.vartrees.wind.wind_ramp_abs:
            wind.append('  wind_ramp_abs %3.6f %3.6f %3.6f %3.6f ;'%
                        (wind_ramp[0],
                         wind_ramp[1],
                         wind_ramp[2],
                         wind_ramp[3]))
        if self.vartrees.wind.scale_time_start > 0:
            wind.append('  scale_time_start        %3.6f ;' %
                        self.vartrees.wind.scale_time_start)
        if self.vartrees.wind.wind_ramp_t1 > 0:
            wind.append('  wind_ramp_factor   %3.6f %3.6f %3.6f %3.6f ;'%
                        (self.vartrees.wind.wind_ramp_t0,
                         self.vartrees.wind.wind_ramp_t1,
                         self.vartrees.wind.wind_ramp_factor0,
                         self.vartrees.wind.wind_ramp_factor1))
        if self.vartrees.wind.iec_gust:
            wind.append('    iec_gust %s %3.6f %3.6f %3.6f %3.6f ;' %
                        (self.vartrees.wind.iec_gust_type,
                         self.vartrees.wind.G_A,
                         self.vartrees.wind_G_phi0,
                         self.vartrees.wind.G_t0,
                         self.vartrees.Wind.G_T))
        if self.vartrees.wind.turb_format == 1:
            wind.append('  begin mann;')
            if self.vartrees.wind.mann.create_turb:
                wind.append('    create_turb_parameters %3.6f %3.6f %3.6f %3.6f %3.6f ;' %
                                  (self.vartrees.wind.mann.L,
                                   self.vartrees.wind.mann.alfaeps,
                                   self.vartrees.wind.mann.gamma,
                                   self.vartrees.wind.mann.seed,
                                   self.vartrees.wind.mann.highfrq_compensation))
            wind.append('    filename_u   %s ;' % (os.path.join(self.turb_directory, self.vartrees.wind.mann.turb_base_name+'_u.bin')))
            wind.append('    filename_v   %s ;' % (os.path.join(self.turb_directory, self.vartrees.wind.mann.turb_base_name+'_v.bin')))
            wind.append('    filename_w   %s ;' % (os.path.join(self.turb_directory, self.vartrees.wind.mann.turb_base_name+'_w.bin')))
            wind.append('    box_dim_u    %i %10.3f ;' % (self.vartrees.wind.mann.box_nu, self.vartrees.wind.mann.box_du))
            wind.append('    box_dim_v    %i %10.3f ;' % (self.vartrees.wind.mann.box_nv, self.vartrees.wind.mann.box_dv))
            wind.append('    box_dim_w    %i %10.3f ;' % (self.vartrees.wind.mann.box_nw, self.vartrees.wind.mann.box_dw))
            wind.append('    std_scaling  %3.6f %3.6f %3.6f ;' % (self.vartrees.wind.mann.std_scaling[0],
                                                                  self.vartrees.wind.mann.std_scaling[1],
                                                                  self.vartrees.wind.mann.std_scaling[2]))
            wind.append('  end mann;')

        if self.vartrees.wind.tower_shadow_method > 0:
            wind.extend(self.write_tower_potential())

        wind.append('end wind;')

        self.master.extend(wind)

    def write_tower_potential(self):

        if hasattr(self.vartrees.wind, 'tower_potential'):
            tower_pot = []
            tower_pot.append('  begin tower_shadow_potential_2 ;')
            tower_pot.append('    tower_mbdy_link %s ;' % self.vartrees.wind.tower_potential.tower_mbdy_link)
            tower_pot.append('    nsec %d; R [m] width [m] (z=rel. distance from node 1 to sec) ;' % self.vartrees.wind.tower_potential.nsec)
            for sec in self.vartrees.wind.tower_potential.sections:
                tower_pot.append('    radius %3.6f %3.6f ;' % (sec[0], sec[1]))

            tower_pot.append('  end tower_shadow_potential_2 ;')
        return tower_pot

    def write_aero(self):

        aero = []
        aero.append('begin aero ;')
        aero.append('  nblades  %d ;' % self.vartrees.aero.nblades)
        aero.append('  hub_vec %s -3 ;' % 'shaft')
        for link in self.vartrees.aero.links:
            aero.append('  link %i mbdy_c2_def %s ;' % (link[0], link[2]))
        aero.append('  ae_filename        ./%s/%s_ae.dat ;' % (self.data_directory, self.case_id))
        aero.append('  pc_filename        ./%s/%s_pc.dat ;' % (self.data_directory, self.case_id))
        aero.append('  induction_method   %i ; 0=none, 1=normal' % self.vartrees.aero.induction_method)
        aero.append('  aerocalc_method    %i ; 0=with aerodynamic, 1=without aerodynamic' % self.vartrees.aero.aerocalc_method)
        aero.append('  aerosections       %i ;' % self.vartrees.aero.aerosections)
        aero.append('  ae_sets            %s ;' % ' '.join(map(str, self.vartrees.aero.ae_sets)))
        aero.append('  tiploss_method     %i ; 0=none, 1=prandtl' % self.vartrees.aero.tiploss_method)
        aero.append('  dynstall_method    %i ; 0=none, 1=stig oeye method,2=mhh method,3=ATEFlap' % self.vartrees.aero.dynstall_method)
        aero.append('    begin dynstall_ateflap ; ')                                                                                                                                          ### tkba: addition for flaps ###
        aero.append('     Ais   %1.2f   %1.2f   %1.2f ;' % (self.vartrees.aero.atef_Ais[0], self.vartrees.aero.atef_Ais[1], self.vartrees.aero.atef_Ais[2]))                                  ### tkba: addition for flaps ###
        aero.append('     Bis   %1.4f   %1.2f   %1.2f ;' % (self.vartrees.aero.atef_Bis[0], self.vartrees.aero.atef_Bis[1], self.vartrees.aero.atef_Bis[2]))                                  ### tkba: addition for flaps ###
        aero.append('     flap   %2.4f   %1.4f   ./%s/%s.ds ; Flap Sec: 1' % (self.vartrees.aero.flap_in, self.vartrees.aero.flap_out, self.data_directory, self.vartrees.aero.ds_filename))  ### tkba: addition for flaps ###
        aero.append('    end dynstall_ateflap ; ')                                                                                                                                            ### tkba: addition for flaps ###
        aero.append('end aero ;')

        self.master.extend(aero)



    def write_aerodrag(self):

        if len(self.vartrees.aerodrag.elements) > 0:
            self.aerodrag.insert(0, 'begin aerodrag ;')
            for i, e in enumerate(self.vartrees.aerodrag.elements):
                self.aerodrag.append('  begin aerodrag_element ;')
                self.aerodrag.append('    mbdy_name %s;' % e.mbdy_name)
                self.aerodrag.append('    aerodrag_sections uniform %i ; distribution of aerodrag calculation points from sec 1 to nsec ;' % e.nsec)
                self.aerodrag.append('    nsec %d; R [m]  Cd [-] width [m] (z=rel. distance from node 1 to sec)' % e.nsec)
                for j in range(e.nsec):
                    self.aerodrag.append('      sec %3.6f %3.6f %3.6f ;' % (e.sections[j][0], e.sections[j][1], e.sections[j][2]))
                self.aerodrag.append('  end aerodrag_element ;')
            self.aerodrag.append('end aerodrag ;')

    def write_structure_res(self):

        structure_out = []
        if self.vartrees.sim.eig_out:
            structure_out.append('beam_output_file_name ./info/%s_beam.dat ;' % self.case_id)
            structure_out.append('body_output_file_name ./info/%s_body.dat ;' % self.case_id)
            structure_out.append('struct_inertia_output_file_name ./info/%s_inertia.dat ;' % self.case_id)
            structure_out.append('body_eigenanalysis_file_name ./info/%s_body_eigs.dat ;' % self.case_id)
            structure_out.append('constraint_output_file_name ./info/%s_constraints.dat ;' % self.case_id)
            structure_out.append('structure_eigenanalysis_file_name ./info/%s_struct_eigs.dat 1 ;' % self.case_id)

        self.structure.extend(structure_out)

    def write_main_bodies(self):
        """
        write all main bodies
        """
        main_bodies = []
        for name in self.vartrees.body_order:

            body = self.vartrees.main_bodies.get_main_body(name)

            main_bodies.append('begin main_body ;')
            if body.copy_main_body is not '':
                main_bodies.append('  name           %s ;' % body.body_name)
                main_bodies.append('  copy_main_body %s ;' % body.copy_main_body)
            else:
                main_bodies.append('  name        %s ;' % body.body_name)
                main_bodies.append('  type        timoschenko ;')
                main_bodies.append('  nbodies     %d ;' % body.nbodies)
                main_bodies.append('  node_distribution     c2_def ;')
                
                if body.damping_type is 'ani':
                    main_bodies.append('  damping_aniso     %s ;' % ' '.join(map(str, body.damping_aniso)))
                else:
                    main_bodies.append('  damping_posdef     %s ;' % ' '.join(map(str, body.damping_posdef)))

                for i in range(len(body.concentrated_mass)):
                    main_bodies.append('  concentrated_mass %s ;' % ' '.join(map(str, body.concentrated_mass[i])))
                main_bodies.append('  begin timoschenko_input;')
                tmpname = ''.join([i for i in body.body_name if not i.isdigit()])
                main_bodies.append('    filename %s ;' % (os.path.join(self.data_directory, self.case_id +'_'+ tmpname + '_st.dat')))
                if body.st_input_type is not 0:
                    main_bodies.append('    becas %d ;' % body.st_input_type)
                main_bodies.append('    set %d %d ;' % (body.body_set[0], body.body_set[1]))
                main_bodies.append('  end timoschenko_input;')
                main_bodies.append('  begin c2_def;')
                main_bodies.append('    nsec %i;' % body.c12axis.shape[0])
                for i in range(body.c12axis.shape[0]):
                    main_bodies.append('    sec %i  %6.6f %6.6f %6.6f %6.6f;' %
                                       (i + 1,
                                        body.c12axis[i, 0],
                                        body.c12axis[i, 1],
                                        body.c12axis[i, 2],
                                        body.c12axis[i, 3]))
                main_bodies.append('  end c2_def;')
                if len(body.beam_structure) > 0:
                    self.write_stfile(body)

            main_bodies.append('end main_body;')

        self.structure.extend(main_bodies)

    def write_orientations(self):

        orientations = []
        orientations.append('begin orientation ;')
        for name in self.vartrees.body_order:
            body = self.vartrees.main_bodies.get_main_body(name)
            for orien in body.orientations:
                if isinstance(orien, HAWC2OrientationBase):
                    orientations.append('  begin base ;')
                    orientations.append('    body %s ;' % body.body_name)
                    orientations.append('    inipos %6.6f %6.6f %6.6f ;' %
                                        (orien.inipos[0],
                                         orien.inipos[1],
                                         orien.inipos[2]))
                    for eulerang in orien.body_eulerang:
                        orientations.append('    body_eulerang %6.6f %6.6f %6.6f ;' %
                                            (eulerang[0],
                                             eulerang[1],
                                             eulerang[2]))
                    orientations.append('  end base ;')

                elif isinstance(orien, HAWC2OrientationRelative):
                    orientations.append('  begin relative ;')
                    fmt = '    body1 '+_get_fmt(orien.body1)+ ';'
                    orientations.append(fmt % (orien.body1[0], orien.body1[1]))
                    fmt = '    body2 '+_get_fmt(orien.body2)+ ';'
                    orientations.append(fmt % (orien.body2[0], orien.body2[1]))
                    for eulerang in orien.body2_eulerang:
                        orientations.append('    body2_eulerang %6.6f %6.6f %6.6f ;' %
                                            (eulerang[0],
                                             eulerang[1],
                                             eulerang[2]))
                    if orien.mbdy2_ini_rotvec_d1[3] != 0.:
                        orientations.append('    body2_ini_rotvec_d1 %3.6f %3.6f %3.6f %3.6f ;' %
                                            (orien.mbdy2_ini_rotvec_d1[0],
                                             orien.mbdy2_ini_rotvec_d1[1],
                                             orien.mbdy2_ini_rotvec_d1[2],
                                             orien.mbdy2_ini_rotvec_d1[3]))
                    orientations.append('  end relative ;')
        orientations.append('end orientation ;')

        self.structure.extend(orientations)

    def write_constraints(self):

        constraints = []

        constraints.append('begin constraint ;')
        for name in self.vartrees.body_order:
            body = self.vartrees.main_bodies.get_main_body(name)
            for con in body.constraints:
                if con.con_type in ['fix0', 'fix2', 'fix3']:
                    constraints.append('  begin %s ;' % con.con_type)
                    constraints.append('    body %s ;' % con.mbdy)
                    if con.disable_at > 0.:
                        constraints.append('    disable_at %s ;' % con.disable_at)
                    if con.con_type in ['fix2', 'fix3']:
                        constraints.append('    dof %i %i %i ;' % (con.dof[0], con.dof[1], con.dof[2]))
                    constraints.append('  end %s ;' % con.con_type)

                elif con.con_type in ['fix1', 'fix4']:
                    constraints.append('  begin %s ;' % con.con_type)
                    fmt = '    body1 '+_get_fmt(con.mbdy1)+ ';'
                    constraints.append(fmt % (con.mbdy1[0], con.mbdy1[1]))
                    fmt = '    body2 '+_get_fmt(con.mbdy2)+ ';'
                    constraints.append(fmt % (con.mbdy2[0], con.mbdy2[1]))
                    constraints.append('  end %s ;' % con.con_type)
                elif 'bearing' in con.con_type:
                    constraints.append('  begin %s ;' % con.con_type)
                    constraints.append('    name %s ; ' % con.name)
                    fmt = '    body1 '+_get_fmt(con.mbdy1)+ ';'
                    constraints.append(fmt % (con.mbdy1[0], con.mbdy1[1]))
                    fmt = '    body2 '+_get_fmt(con.mbdy2)+ ';'
                    constraints.append(fmt % (con.mbdy2[0], con.mbdy2[1]))
                    constraints.append('    bearing_vector %i %3.6f %3.6f %3.6f ;' %
                                       (con.bearing_vector[0],
                                        con.bearing_vector[1],
                                        con.bearing_vector[2],
                                        con.bearing_vector[3]))
                    if con.con_type == 'bearing3':
                        constraints.append('    omegas %3.6f %3.6f %3.6f ;' % (con.omegas))
                    else:
                        if con.disable_at > 0:
                            constraints.append('    disable_at %s ;' % con.disable_at)
                    constraints.append('  end %s ;' % con.con_type)

        constraints.append('end constraint ;')

        self.structure.extend(constraints)

    def write_control_dll(self):
        """
        """

        # Still missing some parameters for the DTUBasicControllerVT

        self.controlinp.append('begin dll ;')
        self.controlinp.append('  begin type2_dll;')
        self.controlinp.append('    name risoe_controller ;')
        self.controlinp.append('    filename  %s/risoe_controller.dll ;'%self.control_directory)
        if self.vartrees.dlls.risoe_controller.dll_init.rotorspeed_gs == 1:
            self.controlinp.append('    dll_subroutine_init init_regulation_2 ;')
        else:
            self.controlinp.append('    dll_subroutine_init init_regulation ;')
        self.controlinp.append('    dll_subroutine_update update_regulation ;')
        self.controlinp.append('    arraysizes_init  50 1 ;')
        self.controlinp.append('    arraysizes_update  9 100 ;')
        self.controlinp.append('    begin init ;')
        self.controlinp.append(';     Overall parameters')
        self.controlinp.append('      constant   1 %3.6e ; Rated power [kW] ' %
                               (self.vartrees.dlls.risoe_controller.dll_init.
                               ratedPower * 1e-3))
        self.controlinp.append('      constant   2 %3.3e ; Minimum rotor '\
                               'speed [rad/s]' %
                               (self.vartrees.dlls.risoe_controller.dll_init.minRPM * 2*np.pi/60))
        self.controlinp.append('      constant   3 %3.6e ; Rated rotor speed [rad/s]' %
                               (self.vartrees.dlls.risoe_controller.dll_init.maxRPM* 2*np.pi/60))
        self.controlinp.append('      constant   4 %3.6e ; Maximum allowable generator torque [Nm]' %
                                self.vartrees.dlls.risoe_controller.dll_init.maxTorque)
        self.controlinp.append('      constant   5 %3.6e ; Minimum pitch angle, theta_min [deg]'  %
                                self.vartrees.dlls.risoe_controller.dll_init.minPitch)
        self.controlinp.append('                         ; if |theta_min|>90, then a table of <wsp,theta_min> is read ;')
        self.controlinp.append('                         ; from a file named "wptable.n", where n=int(theta_min)')
        self.controlinp.append('      constant   6  %3.6e ; Maximum pitch angle [deg]' %
                               self.vartrees.dlls.risoe_controller.dll_init.maxPitch)
        self.controlinp.append('      constant   7  %3.6e ; Maximum pitch velocity operation [deg/s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.maxPitchSpeed)
        self.controlinp.append('      constant   8  %3.6e ; Frequency of generator speed filter [Hz]'   %
                               self.vartrees.dlls.risoe_controller.dll_init.generatorFreq)
        self.controlinp.append('      constant   9  %3.6e ; Damping ratio of speed filter [-]' %
                               self.vartrees.dlls.risoe_controller.dll_init.generatorDamping)
        self.controlinp.append('      constant  10  %3.6e ; Frequency of free-free DT torsion mode [Hz], if zero no notch filter used' %
                               self.vartrees.dlls.risoe_controller.dll_init.ffFreq)
        self.controlinp.append(';     Partial load control parameters')
        self.controlinp.append('      constant  11   %3.6e ; Optimal Cp tracking K factor [kNm/(rad/s)^2];' %
                               self.vartrees.dlls.risoe_controller.dll_init.Qg)
        self.controlinp.append('                           ; Qg=K*Omega^2, K=eta*0.5*rho*A*Cp_opt*R^3/lambda_opt^3')
        self.controlinp.append('      constant  12   %3.6e ; Proportional gain of torque controller [Nm/(rad/s)]' %
                               self.vartrees.dlls.risoe_controller.dll_init.pgTorque)
        self.controlinp.append('      constant  13   %3.6e ; Integral gain of torque controller [Nm/rad] ' %
                               self.vartrees.dlls.risoe_controller.dll_init.igTorque)
        self.controlinp.append('      constant  14   %3.6e ; Differential gain of torque controller [Nm/(rad/s^2)]' %
                               self.vartrees.dlls.risoe_controller.dll_init.dgTorque)
        self.controlinp.append(';     Full load control parameters')
        self.controlinp.append('      constant  15   %i ; Generator control switch [1=constant power, 2=constant torque]' %
                               self.vartrees.dlls.risoe_controller.dll_init.generatorSwitch)
        self.controlinp.append('      constant  16   %3.6e  ; Proportional gain of pitch controller [rad/(rad/s)]' %
                               self.vartrees.dlls.risoe_controller.dll_init.pgPitch)
        self.controlinp.append('      constant  17   %3.6e  ; Integral gain of pitch controller [rad/rad]' %
                               self.vartrees.dlls.risoe_controller.dll_init.igPitch)
        self.controlinp.append('      constant  18   %3.6e    ; Differential gain of pitch controller [rad/(rad/s^2)]' %
                               self.vartrees.dlls.risoe_controller.dll_init.dgPitch)
        self.controlinp.append('      constant  19   %3.6e ; Proportional power error gain [rad/W]' %
                               self.vartrees.dlls.risoe_controller.dll_init.prPowerGain)
        self.controlinp.append('      constant  20   %3.6e ; Integral power error gain [rad/(Ws)]'  %
                               self.vartrees.dlls.risoe_controller.dll_init.intPowerGain)
        self.controlinp.append('      constant  21   %3.6e ; Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]' %
                               self.vartrees.dlls.risoe_controller.dll_init.KK1)
        self.controlinp.append('      constant  22   %3.6e ; Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2] &' %
                               self.vartrees.dlls.risoe_controller.dll_init.KK2)
        self.controlinp.append('                           ; (if zero, KK1 = pitch angle at double gain)')
        self.controlinp.append('      constant  23   %3.6e ; Relative speed for double nonlinear gain [-]' %
                               self.vartrees.dlls.risoe_controller.dll_init.nlGainSpeed)
        self.controlinp.append(';     Cut-in simulation parameters')
        self.controlinp.append('      constant  24   %3.6e ; Cut-in time [s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.cutin_t0)
        self.controlinp.append('      constant  25   %3.6e ; Time delay for soft start of torque [1/1P]' %
                               self.vartrees.dlls.risoe_controller.dll_init.softDelay)
        self.controlinp.append(';     Cut-out simulation parameters')
        self.controlinp.append('      constant  26   %3.6e ; Cut-out time [s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.stop_t0)
        self.controlinp.append('      constant  27   %3.6e ; Time constant for 1st order filter lag of torque cut-out [s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.TorqCutOff)
        self.controlinp.append('      constant  28   1     ; Stop type [1=linear two pitch speed stop, 2=exponential pitch speed stop]')
        self.controlinp.append('      constant  29   %3.6e ; Time delay for pitch stop 1 [s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.PitchDelay1)
        self.controlinp.append('      constant  30   %3.6e ; Maximum pitch velocity during stop 1 [deg/s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.PitchVel1)
        self.controlinp.append('      constant  31   %3.6e ; Time delay for pitch stop 2 [s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.PitchDelay2)
        self.controlinp.append('      constant  32   %3.6e ; Maximum pitch velocity during stop 2 [deg/s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.PitchVel2)
        self.controlinp.append(';     Expert parameters (keep default values unless otherwise given)')
        self.controlinp.append('      constant  33   0.5 ; Lower angle above '
                               'lowest minimum pitch angle for switch [deg]')
        self.controlinp.append('      constant  34   0.5 ; Upper angle above '
                               ' lowest minimum pitch angle for switch [deg], '
                               ' if equal then hard switch')
        self.controlinp.append('      constant  35  95.0 ; Ratio between '
                               ' filtered speed and reference speed for '
                               'fully open torque limits [percantage]')
        self.controlinp.append('      constant  36   5.0 ; Time constant of '
                               '1st order filter on wind speed used for '
                               'minimum pitch [1/1P]')
        self.controlinp.append('      constant  37   5.0 ; Time constant of '
                               '1st order filter on pitch angle used for gain '
                               'scheduling [1/1P]')
        self.controlinp.append(';     Drivetrain damper')
        self.controlinp.append('      constant  38   0.0 ; Proportional gain '
                               'of active DT damper [Nm/(rad/s)], requires '
                               'frequency in input 10')
        self.controlinp.append(';     Over speed')
        self.controlinp.append('      constant  39  %3.6e ; Over speed '
                               'percentage before initiating shut-down' %
                               self.vartrees.dlls.risoe_controller.dll_init.overspeed_limit)
        self.controlinp.append(';    Additional non-linear pitch control term')
        self.controlinp.append('      constant  40  0.0 ; Err0 [rad/s] ')
        self.controlinp.append('      constant  41  0.0 ; ErrDot0 [rad/s^2]')
        self.controlinp.append('      constant  42  0.0 ; PitNonLin1 [rad/s]')

        if self.vartrees.dlls.risoe_controller.dll_init.rotorspeed_gs == 1:
            self.controlinp.append('      constant 43 0.0 ;')
            self.controlinp.append('      constant 44 0.0 ;')
            self.controlinp.append('      constant 45 0.0 ;')
            self.controlinp.append('      constant 46 %f ; kp_speed' %
                                   self.vartrees.dlls.risoe_controller.dll_init.Kp2)
            self.controlinp.append('      constant 47 %f ; invkk1_speed' %
                                   self.vartrees.dlls.risoe_controller.dll_init.Ko1)
            self.controlinp.append('      constant 48 %f ; invkk2_speed' %
                                   self.vartrees.dlls.risoe_controller.dll_init.Ko2)

        self.controlinp.append('   end init ;')
        self.controlinp.append(';')
        self.controlinp.append('    begin output ;')
        self.controlinp.append('      general time ; [s] ')
        self.controlinp.append('      constraint bearing1 shaft_rot 1 only 2 ; Drivetrain speed [rad/s]')
        for i in range(self.vartrees.aero.nblades):
            self.controlinp.append('      constraint bearing2 pitch%d 1 only 1; [rad]' %
                                   (i + 1))
        self.controlinp.append('      wind free_wind 1 0.0 0.0 %3.6f ; global coords at hub height' %
                               self.vartrees.wind.center_pos0[2])
        self.controlinp.append('     general constant 0.0      ; Pitch rate from external system [rad/s]')
        self.controlinp.append('    end output;')
        self.controlinp.append('  end type2_dll;')
        self.controlinp.append(';')
        self.vartrees.dlls.risoe_controller.dll_init.active = True
        if self.vartrees.dlls.risoe_controller.dll_init.active:
            self.controlinp.append('  begin type2_dll;')
            self.controlinp.append('    name generator_servo ;')
            self.controlinp.append('    filename  %s/generator_servo.dll ;'%self.control_directory)
            self.controlinp.append('    dll_subroutine_init init_generator_servo ;')
            self.controlinp.append('    dll_subroutine_update update_generator_servo ;')
            self.controlinp.append('    arraysizes_init  6 1 ;')
            self.controlinp.append('    arraysizes_update  3 6 ;')
            self.controlinp.append('    begin init ;')
            self.controlinp.append('      constant 1  20.0  ; Frequency of genertor 2nd order control model [Hz] ')
            self.controlinp.append('      constant 2  0.9 ; Damping ratio of genertor 2nd order control model [-]')
            self.controlinp.append('      constant 3 %3.6e ; Maximum allowable LSS torque (pull-out torque) [Nm]' %
                                   self.vartrees.dlls.risoe_controller.dll_init.maxTorque)
            self.controlinp.append('      constant 4 %3.6e  ; Generator efficiency [-]' %
                                   self.vartrees.dlls.risoe_controller.dll_init.generatorEfficiency)
            self.controlinp.append('      constant 5 1.0 ; Gear ratio [-]')
            self.controlinp.append('    end init ;')
            self.controlinp.append(';')
            self.controlinp.append('    begin output;')
            self.controlinp.append('      general time                          ;   Time [s]')
            self.controlinp.append('      dll inpvec 1 1                        ;   Electrical torque reference [Nm]')
            self.controlinp.append('      constraint bearing1 shaft_rot 1 only 2;   Generator LSS speed [rad/s]')
            self.controlinp.append('    end output;')
            self.controlinp.append(';')
            self.controlinp.append('    begin actions;')
            self.controlinp.append('      mbdy moment_int shaft 1 -3 shaft towertop 2 ;   Generator LSS torque [Nm]')
            self.controlinp.append('    end actions;')
            self.controlinp.append('  end type2_dll;')
        self.controlinp.append(';')
        self.controlinp.append('  begin type2_dll;')
        self.controlinp.append('    name servo_with_limits ;')
        self.controlinp.append('    filename  %s/servo_with_limits.dll ;'%self.control_directory)
        self.controlinp.append('    dll_subroutine_init init_servo_with_limits ;')
        self.controlinp.append('    dll_subroutine_update update_servo_with_limits ;')
        self.controlinp.append('    arraysizes_init  7 1 ;')
        self.controlinp.append('    arraysizes_update  4 9 ;')
        self.controlinp.append('   begin init ;')
        self.controlinp.append('      constant 1  3     ; 1: Number of blades [-]')
        self.controlinp.append('      constant 2  1.0   ; 2: Filter frequency [Hz]')
        self.controlinp.append('      constant 3  0.7   ; 3: Filter damping ratio [-]')
        self.controlinp.append('      constant 4  %3.6f   ; 4: Max. pitch speed [deg/s]' %
                               self.vartrees.dlls.risoe_controller.dll_init.maxServoPitchSpeed)
        self.controlinp.append('      constant 5  %3.6f   ; 5: Max. pitch acceleration [deg/s^2]' %
                               self.vartrees.dlls.risoe_controller.dll_init.maxServoPitchAcc)
        self.controlinp.append('      constant 6  %3.6f   ; 6: Min. pitch angle [deg] ' %
                               self.vartrees.dlls.risoe_controller.dll_init.minServoPitch)
        self.controlinp.append('      constant 7 90.0   ; 7: Max. pitch angle [deg]')
        self.controlinp.append('   end init ;')
        self.controlinp.append('    begin output;')
        self.controlinp.append('      general time       ;  1: Time                         [s]')
        if self.vartrees.dlls.risoe_controller.dll_init.active:
            for i in range(self.vartrees.rotor.nblades):
                self.controlinp.append('     dll inpvec 1 %d     ;  %d: Pitch%d demand angle          [rad]' % (i+2,i+2,i+1))

        if self.vartrees.dlls.risoe_controller.dll_init.FixedPitch:
            for i in range(self.vartrees.rotor.nblades):
                self.controlinp.append('    general step 3.0 0.0  [PitchAngle%d] ;' % (i +1 ))

        self.controlinp.append('    end output;')
        self.controlinp.append(';')
        self.controlinp.append('    begin actions;')
        for i in range(self.vartrees.rotor.nblades):
            self.controlinp.append('      constraint bearing2 angle pitch%d ; Angle pitch%d bearing    [rad]' % (i+1,i+1))
        self.controlinp.append('    end actions;')
        self.controlinp.append('  end type2_dll;')
        self.controlinp.append(';')
        self.controlinp.append('end dll;')

    def update_c12axis(self):

        self.vartrees.main_bodies.blade1.c12axis = self.vartrees.blade_ae.c12axis.copy()
        self.vartrees.main_bodies.blade1.beam_structure = self.vartrees.blade_structure

    def write_output_fromfile(self):

        self.sensors.append('begin output ;')
        self.sensors.append(' filename %s ;' %
                            (os.path.join(self.res_directory,self.case_id)))
        self.sensors.append(' time %3.6f %3.6f ;' %
                            (self.vartrees.output.time_start,
                             self.vartrees.sim.time_stop))
        self.sensors.append(' data_format %s ;' %
                             self.vartrees.output.out_format)
        self.sensors.append(' buffer 1 ;')

        for i in range(len(self.vartrees.output.sensor_list)):
            self.sensors.append(' ' + self.vartrees.output.sensor_list[i] + ' ;')
        self.sensors.append('end output ;')

    def write_output(self):

        # THIS IS TEMPORARY
        # WE NEED TO MAKE METHODS FOR USERS TO EASILY ADD SENSORS

        self.sensors.append('begin output ;')
        self.sensors.append(' filename %s ;' % (os.path.join(self.res_directory, self.case_id)))
        self.sensors.append(' time %3.6f %3.6f ;' % (self.vartrees.output.time_start, self.vartrees.sim.time_stop))
        self.sensors.append(' data_format %s ;'    % self.vartrees.output.out_format)
        self.sensors.append(' buffer 1 ;')
        self.sensors.append(' general time ;')
        self.sensors.append(' constraint bearing1 shaft_rot 2 ; angle and angle velocity')
        self.sensors.append(' constraint bearing2 pitch1 5 ;    angle and angle velocity')
        self.sensors.append(' constraint bearing2 pitch2 5 ;    angle and angle velocity')
        self.sensors.append(' constraint bearing2 pitch3 5 ;    angle and angle velocity')
        self.sensors.append(' aero omega ;')
        self.sensors.append(' aero torque ;')
        self.sensors.append(' aero power ;')
        self.sensors.append(' aero thrust ;')
        self.sensors.append(' wind free_wind 1 0.0 0.0 %3.6f ; local wind at fixed position: coo (1=global,2=non-rotation rotor coo.), pos x, pos y, pos z' % -self.rotor.hub_height)
        self.sensors.append(' mbdy momentvec tower 1 1  tower # tower base flange ;')
        self.sensors.append(' mbdy momentvec towertop 1 2 towertop # yaw bearing ;')
        self.sensors.append(' mbdy forcevec  towertop 1 2 towertop # yaw bering ;')
        shaft = self.get_main_body('shaft')
        self.sensors.append(' mbdy momentvec shaft %d 1  shaft # main bearing ;' % (shaft.c12axis.shape[0]-1))
        self.sensors.append(' mbdy momentvec blade1 1  1 blade1 # blade 1 root ;')
        self.sensors.append(' mbdy forcevec blade1 1  1 blade1 # blade 1 root ;')
        self.sensors.append(' mbdy momentvec blade1 10 1 local # blade 1 50 percent local e coo ;')
        self.sensors.append('end output ;')

    def write_sensors(self):

        # OUTDATED - DOESN'T WORK

        path = self.sensor_htcfile
        fid = open(path, 'w')
        bladeid = 0
        for n, blade in self.blade_sensors.iteritems():
            bladeid += 1
            for name, sensor in blade.iteritems():
                self.master.append('begin output;')
                self.master.append('    filename %s;' % (os.path.join(self.res_directory, self.blade_resfile) + '_' + name + '_' + str(bladeid)))
                self.master.append('    time %s %s ;' % (self.output_tstart, self.output_tstop))
                self.master.append('    data_format %s ;' % self.output_format)
                self.master.append('    buffer %s ;' % self.output_buffer)

                if sensor['name'] == 'azimuth':
                    self.master.append('    aero azimuth 1;')
                elif sensor['name'] in sensor_type1:
                    for point in sensor['sections']:
                        self.master.append('    aero %s %i %f ;' % (sensor['name'], bladeid, point))
                elif sensor['name'] in sensor_type2:
                    for point in sensor['sections']:
                        self.master.append('    aero %s %i %i %f ;' % (sensor['name'], bladeid, sensor['dof'], point))
                elif sensor['name'] in sensor_type3:
                    for point in sensor['sections']:
                        self.master.append('    aero %s %i %i %i %f ;' % (sensor['name'], sensor['coordsys'], bladeid, sensor['dof'], point))
                self.master.append('end output;')

        self.master.append('begin output;')
        self.master.append('    filename %s;' % (os.path.join(self.res_directory, self.rotor_resfile)))
        self.master.append('    time %s %s ;' % (self.output_tstart, self.output_tstop))
        self.master.append('    data_format  %s ;' % self.output_format)
        self.master.append('    buffer %s ;' % self.output_buffer)
        # self.master.append('    general time ;')
        i = 0
        for name, sensor in self.rotor_sensors.iteritems():
            self.master.append('    aero %s ;' % sensor['name'])
            self.rotor_sensors[name]['pos'] = i
            i += 1
        self.master.append('end output;')
        fid.close()

    def calculate_c12axis(self):
        """
        compute the 1/2 chord axis based on the blade axis and chordwise rotation point

        nb: this examples only works for straight blades!
        """

        # The HAWC2 blade axis is defined using the 1/2 chord points
        b = self.vartrees.blade_geom
        c12axis = np.zeros((b.main_axis.shape[0], 4))
        for i in range(b.main_axis.shape[0]):
            xc12 = (0.5 - b.p_le[i]) * b.chord[i] * np.cos(b.rot_z[i] * np.pi / 180.)
            yc12 = - (0.5 - b.p_le[i]) * b.chord[i] * np.sin(b.rot_z[i] * np.pi / 180.)
            c12axis[i, 0] = -(b.main_axis[i, 0] + xc12)
            c12axis[i, 1] = b.main_axis[i, 1] + yc12
            c12axis[i, 2] = b.main_axis[i, 2] - b.main_axis[0, 2]
        c12axis[:,3] = b.rot_z
        return c12axis

    def write_aefile(self):

        path = os.path.join(self.data_directory, self.case_id + '_ae.dat')
        write_aefile(path, self.vartrees.blade_ae)

        ###### QUICK FIX: COPY DS FILE FROM TOP DATA FOLDER ###################
        try:
            import shutil        
            self.source_path = os.path.join(self.basedir + '/data', self.vartrees.aero.ds_filename + '.ds')     ### tkba: addition for flaps ###  
            self.destin_path = os.path.join(self.data_directory, self.vartrees.aero.ds_filename + '.ds')        ### tkba: addition for flaps ###          
            shutil.copy (self.source_path, self.destin_path)                                                    ### tkba: addition for flaps ###
        except:
            pass
        #######################################################################

    def write_stfile(self, body):

        tmpname = ''.join([i for i in body.body_name if not i.isdigit()])
        path = os.path.join(self.data_directory, self.case_id + '_' + tmpname + '_st.dat')
        write_stfile(path, body, self.case_id)

    def write_pcfile(self):
        
        path = os.path.join(self.data_directory, self.case_id + '_pc.dat')
        write_pcfile(path, self.vartrees.airfoildata)


class HAWC2AeroInputWriter(HAWC2InputWriter):

    def write_main_bodies(self):
        """
        specific writer of the HAWC2aero blade c2_axis
        """

        pass
        # fid = open(path, 'w')
        # lines = []
        # lines.append('begin blade_c2_def;')
        # lines.append('    blade 0 %i;' % b.main_axis.shape[0])
        # for i in range(b.main_axis.shape[0]):
        #     lines.append('        sec %i  %f %f %f %f ;' % (i + 1,
        #                                                     c12axis[i, 0],
        #                                                     c12axis[i, 1],
        #                                                     c12axis[i, 2],
        #                                                     c12axis[i, 3]))
        # lines.append('end blade_c2_def;')
        # for line in lines:
        #     self.master.append(line+'')
        # fid.close()


class HAWC2SInputWriter(HAWC2InputWriter):

    h2s = List()
    set_tsr_flag = Bool(False, iotype='in', desc='Manually set omega according to TSR')
    include_torsiondeform = Int(0, desc='flag for including elasticity (0 = False, 1 = True')
    wsp_cases = List(iotype='in', desc='Wind speeds for which to run HAWC2S')
    cases = List(iotype='in', desc='Other cases to run')

    def write_all(self):
        self.configure_wt()
        self.write_aero()

        turb_format = self.vartrees.wind.turb_format
        self.vartrees.wind.turb_format = 0
        self.write_wind()
        self.vartrees.wind.turb_format = turb_format

        self.write_structure_res()
        self.update_c12axis()
        self.write_main_bodies()
        self.write_orientations()
        self.write_constraints()
        self.structure.insert(0, 'begin new_htc_structure;')
        self.structure.append('end new_htc_structure;')

        self.master.extend(self.structure)
        self.write_hawcstab2()

    def execute(self):

        self.h2s = []
        self.master = []
        self.structure = []
        self.aerodrag = []
        self.controlinp = []
        self.sensors = []

        if not os.path.exists(self.data_directory):
            os.mkdir('data')
        self.case_idout = self.case_id

        self.write_all()
        self.write_master()
        self.write_pcfile()
        self.write_aefile()

    def write_hawcstab2(self):

        self.write_hawcstab2_structure()
        self.write_operational_data_file()
        self.write_h2s_operational_data()
        self.write_h2s_control()
        self.write_h2s_commands()
        self.h2s.insert(0, 'begin hawcstab2 ;')
        self.h2s.append('end hawcstab2 ;')
        self.master.extend(self.h2s)

    def write_h2s_commands(self):

        for name in self.vartrees.h2s.commands:
            if name == 'compute_optimal_pitch_angle':
                self.h2s.append('compute_optimal_pitch_angle use_operational_data ;')

            elif name == 'compute_steady_states' :
                self.h2s.append('compute_steady_states %s %s %s %s ;' %
                                (self.vartrees.h2s.options.bladedeform,
                                 self.vartrees.h2s.options.tipcorrect,
                                 self.vartrees.h2s.options.induction,
                                 self.vartrees.h2s.options.gradients))

            elif name == 'compute_steadystate' :
                self.h2s.append('compute_steadystate %s %s %s %s ;' %
                                (self.vartrees.h2s.options.bladedeform,
                                 self.vartrees.h2s.options.tipcorrect,
                                 self.vartrees.h2s.options.induction,
                                 self.vartrees.h2s.options.gradients))
            elif name == 'compute_stability_analysis' :
                self.h2s.append('compute_stability_analysis %s %s %i %3.6f'
                                 ' %3.6f %3.6f %3.6f %s ;'
                                % (self.vartrees.h2s.options.matrixwriteout,
                                   self.vartrees.h2s.options.eigenvaluewriteout,
                                   self.vartrees.h2s.options.number_of_modes,
                                   self.vartrees.h2s.options.maximum_damping,
                                   self.vartrees.h2s.options.minimum_frequency,
                                   self.vartrees.h2s.options.zero_pole_threshold,
                                   self.vartrees.h2s.options.aero_deflect_ratio,
                                   self.vartrees.h2s.options.frequencysorting))

            elif name == 'compute_aeroservoelastic':
                self.h2s.append('compute_aeroservoelastic %s %s %i %3.6f %3.6f'
                                ' %3.6f %3.6f %s ;' %
                                (self.vartrees.h2s.options.matrixwriteout,
                                 self.vartrees.h2s.options.eigenvaluewriteout,
                                 self.vartrees.h2s.options.number_of_modes,
                                 self.vartrees.h2s.options.maximum_damping,
                                 self.vartrees.h2s.options.minimum_frequency,
                                 self.vartrees.h2s.options.zero_pole_threshold,
                                 self.vartrees.h2s.options.aero_deflect_ratio,
                                 self.vartrees.h2s.options.frequencysorting))
            elif name == 'compute_controller_input':
                self.h2s.append('compute_controller_input ;')
            elif name == 'save_beam_data':
                self.h2s.append('save_beam_data ;')
            elif name == 'save_blade_geometry':
                self.h2s.append('save_blade_geometry ;')
            elif name == 'save_aero_point_data':
                self.h2s.append('save_aero_point_data ;')
            elif name == 'save_profile_coeffs':
                self.h2s.append('save_profile_coeffs ;')
            elif name == 'save_power':
                self.h2s.append('save_power ;')
            elif name == 'save_induction':
                self.h2s.append('save_induction ;')
            elif name == 'save_ol_matrices':
                self.h2s.append('save_ol_matrices ;')
            elif name == 'save_cl_matrices_all':
                if self.vartrees.h2s.options.vloc_out:
                    self.h2s.append('save_cl_matrices_all vloc_out;')
                else:
                    self.h2s.append('save_cl_matrices_all ;')
            elif name == 'compute_structural_modal_analysis':
                if self.vartrees.h2s.options.blade_only:
                    self.h2s.append('compute_structural_modal_analysis '
                                    'bladeonly %i;' %
                                    self.vartrees.h2s.options.number_of_modes)
                else:
                    self.h2s.append('compute_structural_modal_analysis '
                                    'nobladeonly %i;' %
                                    self.vartrees.h2s.options.number_of_modes)
            elif 'steady_state_convergence_limits' in name:
                self.h2s.append(name+' ;')
            elif 'degrees_of_freedom' in name:
                self.h2s.append(name+' ;')
            elif name == 'print_full_precision':
                self.h2s.append('print_full_precision ;')

            elif name == 'basic_dtu_we_controller':
                init = self.vartrees.dlls.risoe_controller.dll_init
                if init.generatorSwitch == 2:
                    generatorSwitch = 0
                elif init.generatorSwitch == 1:
                    generatorSwitch = 1
                self.h2s.append('basic_dtu_we_controller %1.9e %1.9e %1.9e '
                                '%1.9e %1.9e %1.9e %1.9e %1.9e %1.9e %1.9e '
                                '%i %1.9e %1.9e %1.9e;' %
                                (init.pgTorque, init.igTorque, init.Qg,
                                 init.pgPitch, init.igPitch,
                                 init.KK1, init.KK2,
                                 init.generatorFreq, init.generatorDamping,
                                 init.ffFreq, generatorSwitch,
                                 init.Kp2, init.Ko1, init.Ko2))

    def write_operational_data_file(self):

        h2s = self.vartrees.h2s
        ctrl = self.vartrees.dlls.risoe_controller.dll_init

        # HAWCStab2 will compute the operational point for us
        wsp = []
        pitch = []
        rpm = []
        if 'compute_optimal_pitch_angle'  in self.vartrees.h2s.commands \
            and len(h2s.wsp_cases) > 0:

            ctrl.Vin = h2s.wsp_cases[0]
            ctrl.Vout = h2s.wsp_cases[-1]
            ctrl.nV = len(h2s.wsp_cases)

        # operational point is interpolated from the .opt file
        elif h2s.wsp_curve.shape[0] > 0:

            self._logger.info('Operational data wsp_cases %d'%len(h2s.wsp_cases))
            for w in h2s.wsp_cases:
                if self.set_tsr_flag:
                    minRPM = ctrl.minRPM / ctrl.gearRatio
                    maxRPM = ctrl.maxRPM / ctrl.gearRatio
                    omega = ctrl.designTSR * w / self.vartrees.blade_ae.radius
                    r = max(minRPM, min(maxRPM, omega * 60 / (2. * np.pi)))
                    self._logger.info('Setting RPM according to designTSR. wsp = %f TSR = %f RPM = %f' % (w, ctrl.designTSR, r))
                else:
                    r = np.interp(w, h2s.wsp_curve, h2s.rpm_curve)
                p = np.interp(w, h2s.wsp_curve, h2s.pitch_curve)
                wsp.append(w)
                pitch.append(p)
                rpm.append(r)
                self._logger.info('adding case %f %f %f' % (w, p, r))
        for case in h2s.cases:
            try:
                wsp.append(case['wsp'])
                pitch.append(case['pitch'])
                rpm.append(case['rpm'])
                self._logger.info('adding case %f %f %f' % (case['wsp'], case['pitch'], case['rpm']))

            except:
                raise RuntimeError('wrong inputs in case')

        if len(wsp) > 0:
            data = np.array([wsp, pitch, rpm]).T
            ctrl.Vin = wsp[0]
            ctrl.Vout = wsp[-1]
            ctrl.nV = len(wsp)
        else:
            data = np.array([h2s.wsp_curve, h2s.pitch_curve, h2s.rpm_curve]).T

        fid = open(self.case_id + '.opt', 'w')
        fid.write('%i Wind speed [m/s]          Pitch [deg]     Rot. speed [rpm]\n' % data.shape[0])
        np.savetxt(fid, data)
        fid.close()

    def write_hawcstab2_structure(self):

        if self.from_file:
            self.h2s.append('  begin ground_fixed_substructure ;')
            for name in self.vartrees.h2s.ground_fixed.main_body:
                self.h2s.append('    main_body %s ;' % name)
            if self.vartrees.h2s.ground_fixed.log_decrements[0] != 0:
                self.h2s.append('    log_decrements %3.6f %3.6f ;' %
                                (self.vartrees.h2s.ground_fixed.log_decrements[0],
                                 self.vartrees.h2s.ground_fixed.log_decrements[1]))
            self.h2s.append('  end ground_fixed_substructure ;')

            self.h2s.append('  begin rotating_axissym_substructure ;')
            for name in self.vartrees.h2s.rotating_axissym.main_body:
                self.h2s.append('    main_body %s ;' % name)
            if self.vartrees.h2s.rotating_axissym.log_decrements[0] != 0:
                self.h2s.append('    log_decrements %3.6f %3.6f ;' %
                                (self.vartrees.h2s.rotating_axissym.log_decrements[0],
                                 self.vartrees.h2s.rotating_axissym.log_decrements[1]))
            self.h2s.append('  end rotating_axissym_substructure ;')

            self.h2s.append('  begin rotating_threebladed_substructure ;')
            for name in self.vartrees.h2s.rotating_threebladed.main_body:
                self.h2s.append('    main_body %s ;' % name)
            if self.vartrees.h2s.rotating_threebladed.log_decrements[0] != 0:
                self.h2s.append('    log_decrements %3.6f %3.6f %3.6f %3.6f %3.6f %3.6f ;' %
                                (self.vartrees.h2s.rotating_threebladed.log_decrements[0],
                                 self.vartrees.h2s.rotating_threebladed.log_decrements[1],
                                 self.vartrees.h2s.rotating_threebladed.log_decrements[2],
                                 self.vartrees.h2s.rotating_threebladed.log_decrements[3],
                                 self.vartrees.h2s.rotating_threebladed.log_decrements[4],
                                 self.vartrees.h2s.rotating_threebladed.log_decrements[5]))
            self.h2s.append('    second_order_actuator %s %3.6f %3.6f ;' %
                            (self.vartrees.h2s.second_order_actuator.name,
                             self.vartrees.h2s.second_order_actuator.frequency,
                             self.vartrees.h2s.second_order_actuator.damping))
            self.h2s.append('  end rotating_threebladed_substructure ;')
        else:
            self.h2s.append('  begin ground_fixed_substructure ;')
            self.h2s.append('    main_body tower ;')
            self.h2s.append('    main_body towertop ;')
            self.h2s.append('  end ground_fixed_substructure ;')
            self.h2s.append('  begin rotating_axissym_substructure ;')
            self.h2s.append('    main_body shaft ;')
            self.h2s.append('  end rotating_axissym_substructure ;')
            self.h2s.append('  begin rotating_threebladed_substructure ;')
            self.h2s.append('    main_body hub1 ;')
            self.h2s.append('    main_body blade1 ;')
            self.h2s.append('    second_order_actuator pitch1  100 0.9 ;')
            self.h2s.append('  end rotating_threebladed_substructure ;')

    def write_h2s_control(self):

        self.h2s.append('  begin controller_tuning ;')
        self.h2s.append('    partial_load %3.6f %3.6f; fn [hz], zeta [-]' %
                        (self.vartrees.dlls.risoe_controller.dll_init.poleFreqTorque,
                         self.vartrees.dlls.risoe_controller.dll_init.poleDampTorque))
        self.h2s.append('    full_load %3.6f %3.6f ; fn [hz], zeta [-]' %
                        (self.vartrees.dlls.risoe_controller.dll_init.poleFreqPitch,
                         self.vartrees.dlls.risoe_controller.dll_init.poleDampPitch))
        self.h2s.append('    gain_scheduling %d ; 1 linear, 2 quadratic' %
                        self.vartrees.dlls.risoe_controller.dll_init.gainScheduling)
        self.h2s.append('    constant_power %d ; ' %
                        self.vartrees.dlls.risoe_controller.dll_init.generatorSwitch)
        if len(self.vartrees.h2s.options.regions):
            self.h2s.append('    regions %i %i %i %i ; ' %(
                            self.vartrees.h2s.options.regions[0],
                            self.vartrees.h2s.options.regions[1],
                            self.vartrees.h2s.options.regions[2],
                            self.vartrees.h2s.options.regions[3]))            
        self.h2s.append('  end controller_tuning ;')

        if self.from_file:
            self.h2s.append('  begin controller ;')
            self.h2s.append('    begin input ;')
            for i in range(len(self.vartrees.h2s.ch_list_in.sensor_list)):
                self.h2s.append('      '+self.vartrees.h2s.ch_list_in.sensor_list[i]+';')
            self.h2s.append('    end input ;')
            self.h2s.append('    begin output ;')
            for i in range(len(self.vartrees.h2s.ch_list_out.sensor_list)):
                self.h2s.append('      '+self.vartrees.h2s.ch_list_out.sensor_list[i]+';')
            self.h2s.append('     end output ;')
        else:
            self.h2s.append('  begin controller ;')
            self.h2s.append('    begin input ;')
            self.h2s.append('     constraint bearing1 shaft_rot ;')
            self.h2s.append('     constraint bearing2 pitch1 collective ;')
            self.h2s.append('     constraint bearing2 pitch1 cosine ;')
            self.h2s.append('     constraint bearing2 pitch1 sine ;')
            self.h2s.append('    end input ;')
            self.h2s.append('    begin output ;')
            self.h2s.append('      constraint bearing1 shaft_rot 1 only 2 ; 1')
            self.h2s.append('      constraint bearing2 pitch1 1 only 1 collective ; 2')
            self.h2s.append('      constraint bearing2 pitch1 1 only 1 cosine ; 3')
            self.h2s.append('      constraint bearing2 pitch1 1 only 1 sine ; 4')
            self.h2s.append('      mbdy momentvec hub1 1 0 hub1 only 1 collective ; 5')
            self.h2s.append('      mbdy momentvec hub1 1 0 hub1 only 1 cosine ; 6')
            self.h2s.append('      mbdy momentvec hub1 1 0 hub1 only 1 sine ; 7')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 2 collective ; 17')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 2 cosine ; 18')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 2 sine ; 19')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 1 collective ; 20')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 1 cosine ; 21')
            self.h2s.append('      mbdy momentvec blade1 1 0 blade1 only 1 sine ; 22')
            self.h2s.append('     end output ;')

        self.h2s.append('  end controller ;')

    def write_h2s_operational_data(self):

        self.h2s.append('  operational_data_filename %s ;' %
                        (self.case_id + '.opt'))

        self.h2s.append('  begin operational_data ;')
        self.h2s.append('    windspeed %3.6f %3.6f %d ; cut-in, cut-out, points' %
                        (self.vartrees.dlls.risoe_controller.dll_init.Vin,
                         self.vartrees.dlls.risoe_controller.dll_init.Vout,
                         self.vartrees.dlls.risoe_controller.dll_init.nV))
        self.h2s.append('    genspeed %3.6f %3.6f ;' %
                        (self.vartrees.dlls.risoe_controller.dll_init.minRPM,
                         self.vartrees.dlls.risoe_controller.dll_init.maxRPM))
        self.h2s.append('    gearratio %3.6f ;' %
                        self.vartrees.dlls.risoe_controller.dll_init.gearRatio)
        self.h2s.append('    minpitch %3.6f ;' %
                        self.vartrees.dlls.risoe_controller.dll_init.minPitch)
        self.h2s.append('    opt_lambda %3.6f ;' %
                        self.vartrees.dlls.risoe_controller.dll_init.designTSR)
        self.h2s.append('    maxpow %3.6e ;' %
                        (self.vartrees.dlls.risoe_controller.dll_init.ratedAeroPower))
        self.h2s.append('    prvs_turbine %d ;' %
                        self.vartrees.dlls.risoe_controller.dll_init.prvs_turbine)
        self.h2s.append('    include_torsiondeform %d ;' %
                        self.vartrees.h2s.options.include_torsiondeform)
        self.h2s.append('  end operational_data ;')
