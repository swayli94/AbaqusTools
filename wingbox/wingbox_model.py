'''

'''
import os
import sys
import time
import json
import numpy as np

from AbaqusTools import Model, IS_ABAQUS
from lofting_part import LoftingPart
from rib_part import RibPart
from rib_part import use_rib_face_tie_for_internal_spars
from params import get_parameter_file, load_parameters, get_failure_model

if IS_ABAQUS:
    from abaqus import mdb
    from abaqusConstants import *


ALUMINUM_RIB_TYPES = ('aluminum', 'aluminium', 'aluminum_alloy', 'aluminium_alloy')


class WingboxModel(Model):
    '''
    Class for wingbox model.
    '''
    def __init__(self, name_job, pGeo, pMesh, pRun):
        super(WingboxModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun)
        self.name_job = name_job

    def _get_analysis_type(self):
        '''
        'static', 'buckle', or 'static_buckle', i.e., whether a general static
        step and/or a linear perturbation buckling step is solved.
        '''
        return str(self.pRun.get('analysis_type', 'static')).strip().lower()

    def _get_static_step_name(self):
        if self._get_analysis_type() == 'static_buckle':
            return 'Preload'
        return 'Loading'

    def _get_buckle_step_name(self):
        if self._get_analysis_type() == 'static_buckle':
            return 'Buckling'
        return 'Loading'

    def _get_buckling_base_state(self):
        '''
        Get the base state of the buckling step in a static_buckle analysis.

        'initial': the buckling step is created before the static step.
        'preload': the buckling step is created after the static step.
        '''
        buckling_base_state = str(
            self.pRun.get('buckling_base_state', 'initial')).lower()
        if buckling_base_state not in ('initial', 'preload'):
            raise ValueError(
                'Invalid buckling_base_state: %s' % buckling_base_state)
        return buckling_base_state

    def _get_output_variables(self, name, default):
        '''
        Get user-defined field output variables from `pRun[name]`.

        Parameters
        ----------------
        name: str
            parameter name, i.e., 'output_variables',
            'output_layup_variables', or 'output_buckling_variables'.

        default: tuple of str
            variables used when `pRun[name]` is not defined.

        Notes
        ----------------
        `pRun[name]` can be:

        - a list of variable names, e.g., ["S", "E", "U", "P"]
        - a single variable name, e.g., "U"
        - "PRESELECT" or "ALL", i.e., the Abaqus default or complete set
        - an empty list, i.e., no output request is created
        '''
        if name not in self.pRun:
            return tuple(default)

        variables = self.pRun[name]

        if isinstance(variables, (list, tuple)):
            items = list(variables)
        else:
            items = [variables]

        names = []
        for item in items:
            variable = str(item).strip()
            if variable == '':
                raise ValueError('Empty variable name in pRun.%s.' % name)
            if variable not in names:
                names.append(variable)

        if len(names) == 1 and names[0].upper() in ('PRESELECT', 'ALL'):
            if names[0].upper() == 'PRESELECT':
                return PRESELECT
            return ALL

        return tuple(names)

    def _get_failure_output_variables(self):
        '''
        Field output variables carrying the failure indices of the current
        failure model.  They are evaluated at every section point of the layup,
        and are reduced to one value per element by `postprocess_failure.py`.

        Failure analysis is optional: `pMesh['failure_model']` can be absent,
        null or 'none', and then no failure index is requested at all.
        '''
        failure_model = get_failure_model(self.pMesh)
        if failure_model == 'hashin':
            return ('DMICRT', 'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT')
        if failure_model == 'larc05':
            # The LaRC05 criteria are evaluated by the user subroutine and
            # stored in the user-defined output variables.
            return ('UVARM', )
        return ()

    def _get_default_rib_output_variables(self):
        '''
        Field output variables of the ribs.

        A metallic rib has no composite layup, so it appears in no layup output
        request and its strength can only be checked from its stress field.
        '''
        if str(self.pMesh.get('rib_material_type', 'composite')).lower() in ALUMINUM_RIB_TYPES:
            return ('S', )
        return ()

    def _get_layup_location_options(self):
        '''
        Section points written by the layup field output requests.

        `pRun['layup_output_ply_locations']` can be:

        - 'all' (default): every section point of every ply.
        - 'top_bottom': ply top and ply bottom, i.e., two thirds of the data of
          a three-point Simpson rule.  The failure indices are convex functions
          of a stress state that varies linearly across a ply, so the ply
          extremes are kept unless a stress component changes sign inside the
          ply.
        - 'mid': ply middle only, the smallest and least conservative option.
        '''
        mode = str(self.pRun.get('layup_output_ply_locations', 'all')).lower()
        if mode in ('all', 'all_locations'):
            return {'layupLocationMethod': ALL_LOCATIONS}
        if mode in ('top_bottom', 'top-bottom'):
            return {'layupLocationMethod': SPECIFIED,
                    'outputAtPlyTop': True, 'outputAtPlyMid': False,
                    'outputAtPlyBottom': True}
        if mode == 'mid':
            return {'layupLocationMethod': SPECIFIED,
                    'outputAtPlyTop': False, 'outputAtPlyMid': True,
                    'outputAtPlyBottom': False}
        raise ValueError('Invalid layup_output_ply_locations: %s' % mode)

    def _create_metallic_rib_material(self):
        '''
        Create the isotropic material for metallic ribs.

        `pMesh['rib_material_name']` must be a key of
        `AbaqusTools.materials.MATERIAL_LIBRARY`; add an entry there to use
        another alloy. Units follow the rest of the model, i.e., N, mm, tonne.
        '''
        self.create_material(str(self.pMesh.get('rib_material_name', 'Aluminum-7075')))

    def initialization(self):

        self.model = mdb.models[str(self.name_model)]

        self.create_material_IM785517(elastic_type='ENGINEERING_CONSTANTS')
        if str(self.pMesh.get('rib_material_type', 'composite')).lower() in ALUMINUM_RIB_TYPES:
            self._create_metallic_rib_material()

    def setup_parts(self):
        
        self.lofting = LoftingPart(name_part='lofting',
                        model=self.model, pGeo=self.pGeo, pMesh=self.pMesh)
        self.lofting.build()
        
        self.ribs = []
        for i in range(len(self.pGeo['sections'])):
            rib = RibPart(model=self.model, pGeo=self.pGeo, pMesh=self.pMesh, index_rib=i)
            rib.build()
            self.ribs.append(rib)

    def setup_assembly(self):
        
        a = self.rootAssembly
        p = self.model.parts['lofting']
        a.Instance(name='lofting', part=p, dependent=ON)
        
        for i, rib in enumerate(self.ribs):
            p = self.model.parts[rib.name_part]
            a.Instance(name=rib.name_part, part=p, dependent=ON)
    
    def setup_steps(self):

        analysis_type = self._get_analysis_type()

        buckle_options = {
            'numEigen': self.pRun.get('numEigen', 10),
            'eigensolver': LANCZOS,
        }
        if self.pRun.get('buckling_min_eigen') is not None:
            buckle_options['minEigen'] = float(self.pRun['buckling_min_eigen'])
        if self.pRun.get('buckling_block_size') is not None:
            buckle_options['blockSize'] = int(self.pRun['buckling_block_size'])
        if self.pRun.get('buckling_max_blocks') is not None:
            buckle_options['maxBlocks'] = int(self.pRun['buckling_max_blocks'])

        if analysis_type == 'buckle':
            self.model.BuckleStep(
                name=self._get_buckle_step_name(),
                previous='Initial',
                **buckle_options
                )
            return

        if analysis_type not in ['static', 'static_buckle']:
            raise ValueError('Invalid analysis_type: %s' % analysis_type)

        if self.pRun['nlgeom']:
            nlgeom = ON
        else:
            nlgeom = OFF
        
        buckling_base_state = self._get_buckling_base_state()

        # For production static+buckling analyses, solve the eigenvalue problem
        # directly from the unloaded Initial state.  The resulting eigenvalues
        # are therefore total load multipliers and can be compared directly
        # with the optimisation constraint lambda >= 1.  A subsequent general
        # static step starts from the last general state (Initial); the linear
        # perturbation step does not modify that state.
        if analysis_type == 'static_buckle' and buckling_base_state == 'initial':
            self.model.BuckleStep(
                name=self._get_buckle_step_name(),
                previous='Initial',
                **buckle_options
                )

            static_previous = self._get_buckle_step_name()
        else:
            static_previous = 'Initial'

        self.model.StaticStep(
            name=self._get_static_step_name(),
            previous=static_previous,
            description='Static simulation',
            timePeriod=self.pRun['timePeriod'],
            maxNumInc=self.pRun['maxNumInc'],
            initialInc=self.pRun['initialInc'],
            minInc=self.pRun['minInc'],
            maxInc=self.pRun['maxInc'],
            nlgeom=nlgeom,
            )

        # Retain the former prestressed incremental formulation only as an
        # explicit compatibility option for reproducing archived calculations.
        if analysis_type == 'static_buckle' and buckling_base_state == 'preload':
            self.model.BuckleStep(
                name=self._get_buckle_step_name(),
                previous=self._get_static_step_name(),
                **buckle_options
                )
            
    def setup_interactions(self):

        tie_strategy = str(self.pMesh.get('tie_strategy', 'legacy_sets')).lower()
        if tie_strategy in ('station_nodes', 'station-nodes'):
            self._setup_station_node_interactions()
            return
        if tie_strategy not in ('legacy', 'legacy_sets', 'legacy-sets'):
            raise ValueError('Invalid tie_strategy: %s' % tie_strategy)
        
        n_sections = len(self.lofting.sections)

        #* Tie between covers and ribs
        pairs = []
        for i_section in range(n_sections-1):
            index_ribs = [i_section, i_section+1]
            for side in ['upper', 'lower']:
                nsi = 'lofting'
                nss = 'face_wingbox%d_cover_%s' % (i_section, side)
                for k in index_ribs:
                    nmi = self.ribs[k].name_part
                    nms = 'edge_rib_cover_%s' % (side)
                    pairs.append((nmi, nms, nsi, nss))

        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_cover_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
        #* Tie between spars and ribs
        pairs = []
        for i_section in range(n_sections):
            section = self.lofting.sections[i_section]
            for j in range(section.n_spars):
                nmi = self.ribs[i_section].name_part
                if (
                        use_rib_face_tie_for_internal_spars(self.pMesh)
                        and j > 0
                        and j < section.n_spars - 1):
                    nms = 'face_rib'
                else:
                    nms = 'edge_rib_spar_%d' % (j)
                nsi = 'lofting'
                nss = 'edge_sec%d_spar_%d' % (i_section, j)
                pairs.append((nmi, nms, nsi, nss))
        
        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_spar_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
        #* Tie between stringers and ribs
        pairs = []
        for i_section, section in enumerate(self.lofting.sections):
            nmi = self.ribs[i_section].name_part
            nms = 'face_rib'
            for j in range(section.n_stringers):
                for side in ['upper', 'lower']:
                    for feature in ['web', 'flange']:
                        nsi = 'lofting'
                        nss = 'edge_sec%d_stringer_%d_%s_%s' % (i_section, j, side, feature)
                        pairs.append((nmi, nms, nsi, nss))
                    
        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_stringer_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
    def _create_tie_constraint(self, name_constraint,
                    name_master_instance, name_master_set,
                    name_slave_instance, name_slave_set):
        '''
        Create tie constraint.
        '''
        a = self.rootAssembly
        region1=a.instances[name_master_instance].sets[name_master_set]
        region2=a.instances[name_slave_instance].sets[name_slave_set]
        adjust = ON if self._pmesh_bool('tie_adjust', True) else OFF
        thickness = ON if self._pmesh_bool('tie_thickness', True) else OFF
        self.model.Tie(name=name_constraint, main=region1, secondary=region2,
            positionToleranceMethod=COMPUTED, adjust=adjust,
            tieRotations=ON, thickness=thickness)

    def _pmesh_bool(self, name, default):
        value = self.pMesh.get(name, default)
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized in ('true', 'on', 'yes', '1'):
                return True
            if normalized in ('false', 'off', 'no', '0'):
                return False
            raise ValueError('Invalid boolean value for pMesh.%s: %s' % (name, value))
        return bool(value)

    def _station_interface_set_names(self, i_section, section):
        tag = 'sec%d' % i_section
        names = [
            'edge_%s_cover_upper' % tag,
            'edge_%s_cover_lower' % tag,
        ]
        for j in range(section.n_spars):
            names.append('edge_%s_spar_%d' % (tag, j))
        for j in range(section.n_stringers):
            for side in ('upper', 'lower'):
                for feature in ('web', 'flange'):
                    names.append(
                        'edge_%s_stringer_%d_%s_%s'
                        % (tag, j, side, feature)
                    )
        return names

    def _setup_station_node_interactions(self):
        """Tie each loft station to its rib with one unique node-based secondary."""
        a = self.rootAssembly
        loft_instance = a.instances['lofting']
        owners = {}
        station_tolerance = float(
            self.pMesh.get('tie_station_tolerance', 1.0e-6)
        )
        position_tolerance = float(
            self.pMesh.get('tie_position_tolerance', 0.1)
        )
        if station_tolerance <= 0.0:
            raise ValueError('tie_station_tolerance must be positive.')
        if position_tolerance <= 0.0:
            raise ValueError('tie_position_tolerance must be positive.')

        adjust = ON if self._pmesh_bool('tie_adjust', False) else OFF
        thickness = ON if self._pmesh_bool('tie_thickness', True) else OFF
        tie_rotations = ON if self._pmesh_bool('tie_rotations', True) else OFF
        node_selection = str(
            self.pMesh.get('station_node_selection', 'source_sets')
        ).lower()
        if node_selection not in ('source_sets', 'station_plane'):
            raise ValueError(
                'Invalid station_node_selection: %s' % node_selection)

        for i_section, section in enumerate(self.lofting.sections):
            if len(section.cutouts) != 0:
                raise RuntimeError(
                    'station_nodes tie is not validated for rib cutouts '
                    '(section=%d, cutouts=%d).'
                    % (i_section, len(section.cutouts))
                )

            labels = set()
            candidate_count = 0
            source_set_names = []
            if node_selection == 'station_plane':
                for node in loft_instance.nodes:
                    if abs(
                            float(node.coordinates[2])
                            - float(section.zLE)) <= station_tolerance:
                        labels.add(int(node.label))
                        candidate_count += 1
                source_count = 1
            else:
                source_set_names = self._station_interface_set_names(
                    i_section, section)
                for name_set in source_set_names:
                    if name_set not in loft_instance.sets:
                        raise RuntimeError(
                            'Missing loft station tie set %s at section %d.'
                            % (name_set, i_section)
                        )
                    nodes = loft_instance.sets[name_set].nodes
                    if len(nodes) == 0:
                        raise RuntimeError(
                            'Empty loft station tie set %s at section %d.'
                            % (name_set, i_section)
                        )
                    candidate_count += len(nodes)
                    for node in nodes:
                        if abs(
                                float(node.coordinates[2])
                                - float(section.zLE)) > station_tolerance:
                            raise RuntimeError(
                                'Node %d in %s is off station %d: '
                                'z=%.12g expected %.12g.'
                                % (
                                    int(node.label), name_set, i_section,
                                    float(node.coordinates[2]),
                                    float(section.zLE),
                                )
                            )
                        labels.add(int(node.label))
                source_count = len(source_set_names)

            if not labels:
                raise RuntimeError(
                    'No secondary nodes found for loft station %d.' % i_section
                )

            for label in labels:
                key = ('lofting', label)
                if key in owners:
                    raise RuntimeError(
                        'Loft node %d is secondary at both station %d and station %d.'
                        % (label, owners[key], i_section)
                    )
                owners[key] = i_section

            name_secondary = 'tie_secondary_station_%d' % i_section
            self.create_node_set(
                name_set=name_secondary,
                name_instance='lofting',
                node_labels=tuple(sorted(labels)),
                unsorted=False,
            )

            name_rib = self.ribs[i_section].name_part
            if 'face_rib' not in a.instances[name_rib].surfaces:
                raise RuntimeError(
                    'Missing rib main surface face_rib on %s.' % name_rib
                )
            main_surface = a.instances[name_rib].surfaces['face_rib']
            main_face_count = len(main_surface.faces)
            if main_face_count == 0:
                raise RuntimeError(
                    'Empty rib main surface face_rib on %s.' % name_rib
                )

            self.model.Tie(
                name='tie_station_%d' % i_section,
                main=main_surface,
                secondary=a.sets[name_secondary],
                positionToleranceMethod=SPECIFIED,
                positionTolerance=position_tolerance,
                adjust=adjust,
                tieRotations=tie_rotations,
                constraintEnforcement=NODE_TO_SURFACE,
                thickness=thickness,
            )
            print(
                '>>> STATION_UNION station=%d source_sets=%d raw_memberships=%d '
                'unique_nodes=%d duplicates_removed=%d main_faces=%d'
                % (
                    i_section, source_count, candidate_count,
                    len(labels), candidate_count - len(labels), main_face_count,
                )
            )
    
    def setup_loads(self):
        
        a = self.rootAssembly
        n_sections = len(self.lofting.sections)
        analysis_type = self._get_analysis_type()
        if analysis_type == 'buckle':
            load_step_names = [self._get_buckle_step_name()]
        elif analysis_type == 'static_buckle':
            load_step_names = [self._get_static_step_name(), self._get_buckle_step_name()]
        else:
            load_step_names = [self._get_static_step_name()]
        
        region = a.instances['rib_0'].sets['all']
        self.model.EncastreBC(name='sec0', createStepName='Initial', 
            region=region, localCsys=None)
        
        pressure_start_section = self.pMesh.get('pressure_start_section', 0)
        for load_step_name in load_step_names:
            for side in ['upper', 'lower']:
                for i_section in range(pressure_start_section, n_sections-1):
                    region = a.instances['lofting'].surfaces['face_wingbox%d_cover_%s' % (i_section, side)]
                    self.model.Pressure(name='pressure_%s_%d_%s' % (load_step_name, i_section, side),
                        createStepName=load_step_name, region=region,
                        distributionType=UNIFORM, field='',
                        magnitude=self.pMesh['pressure_%s_cover' % side],
                        amplitude=UNSET)

        # The engine dead weight and the omitted outboard-wing lift are kept as
        # separate named loads for traceability, but intentionally share the
        # existing engine-pylon reference point.  The current reduced model
        # assumes zero eccentricity at that point, so no additional moment is
        # applied for the outboard lift.
        point_loads = []
        if 'engine_weight_force' in self.pMesh:
            point_loads.append(('engine_weight', self.pMesh['engine_weight_force']))
        if 'outboard_lift_force' in self.pMesh:
            application = self.pMesh.get(
                'outboard_lift_application', 'engine_pylon_reference_point')
            if application != 'engine_pylon_reference_point':
                raise ValueError(
                    'Unsupported outboard_lift_application: %s' % application)
            eccentricity = float(self.pMesh.get('outboard_lift_eccentricity_mm', 0.0))
            if abs(eccentricity) > 1.0e-12:
                raise ValueError(
                    'Nonzero outboard lift eccentricity requires an explicit moment model.')
            point_loads.append(('outboard_lift', self.pMesh['outboard_lift_force']))

        if point_loads:
            i_engine = self.pMesh.get('engine_pylon_section', n_sections-1)
            i_spar = self.pMesh.get('engine_pylon_spar_index', 1)
            name_rib = self.ribs[i_engine].name_part
            name_rp = 'RP_engine_pylon'
            point_rp = self.lofting.sections[i_engine].spars[i_spar].get_selection_point(
                feature='spar', side=None)

            self.create_reference_point(point_rp[0], point_rp[1], point_rp[2], name_rp)
            self.create_reference_point_set(name_rp, name_rp)

            self.model.Coupling(name='coupling_engine_pylon_rib',
                controlPoint=a.sets[name_rp],
                surface=a.instances[name_rib].surfaces['face_rib'],
                influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                localCsys=None,
                u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

            for load_step_name in load_step_names:
                for load_name, force_y in point_loads:
                    self.model.ConcentratedForce(
                        name='%s_%s' % (load_name, load_step_name),
                        createStepName=load_step_name,
                        region=a.sets[name_rp],
                        cf2=force_y,
                        amplitude=UNSET,
                        follower=OFF,
                        localCsys=None)
        
    def _create_tip_node_set(self, tolerance=1.0):
        '''
        Assembly-level node set of the wing tip cross-section.

        The wing spans in the global z direction, so the tip is the set of
        nodes with the largest z coordinate (within `tolerance` mm).  The
        history output of this set is what the fast post-processing
        (`extract_results.py`) reads for the tip displacement, instead of
        scanning the full nodal U field of the output database.
        '''
        a = self.rootAssembly

        z_max = None
        for name_instance in a.instances.keys():
            for node in a.instances[name_instance].nodes:
                if z_max is None or node.coordinates[2] > z_max:
                    z_max = node.coordinates[2]
        if z_max is None:
            return None

        labels_of_instance = []
        for name_instance in a.instances.keys():
            labels = [node.label for node in a.instances[name_instance].nodes
                      if node.coordinates[2] >= z_max - tolerance]
            if labels:
                labels_of_instance.append((name_instance, tuple(labels)))

        a.SetFromNodeLabels(name='Tip', nodeLabels=tuple(labels_of_instance))
        print('>>> TIP_NODE_SET %d nodes at z = %.1f mm: %s'
              % (sum([len(item[1]) for item in labels_of_instance]), z_max,
                 {item[0]: len(item[1]) for item in labels_of_instance}))
        return a.sets['Tip']

    def setup_outputs(self):

        analysis_type = self._get_analysis_type()
        output_mode = str(self.pRun.get('output_mode', 'optimisation')).lower()
        if output_mode not in ('full', 'optimisation'):
            raise ValueError('Invalid output_mode: %s' % output_mode)

        #* Default variables of the two output modes
        if output_mode == 'full':
            default_variables = ('S', 'E', 'U', 'P')
            default_layup_variables = ('S', 'E')
        else:
            default_variables = ('U',)
            default_layup_variables = ()

        # Failure criteria are always requested, unless the layup variables are
        # explicitly defined by the user.
        if 'output_layup_variables' not in self.pRun:
            default_layup_variables += self._get_failure_output_variables()

        #* User-defined variables (pRun), which override the output mode
        variables = self._get_output_variables(
            'output_variables', default_variables)
        layup_variables = self._get_output_variables(
            'output_layup_variables', default_layup_variables)
        buckling_variables = self._get_output_variables(
            'output_buckling_variables', ('U',))
        rib_variables = self._get_output_variables(
            'output_rib_variables', self._get_default_rib_output_variables())

        if not variables:
            raise ValueError('pRun.output_variables must not be empty.')
        if not buckling_variables:
            raise ValueError('pRun.output_buckling_variables must not be empty.')

        print('>>> OUTPUT_VARIABLES mode=%s global=%s layup=%s buckling=%s rib=%s'
              % (output_mode, variables, layup_variables, buckling_variables,
                 rib_variables))

        # Abaqus creates PRESELECT field/history requests automatically when a
        # step is created.  If retained, those requests write every converged
        # static increment even when our explicit requests use LAST_INCREMENT.
        # Remove them first so the ODB contains only the controlled outputs
        # below: the first buckling eigenvector and the final static frame.
        for request_name in list(self.model.fieldOutputRequests.keys()):
            del self.model.fieldOutputRequests[request_name]
        for request_name in list(self.model.historyOutputRequests.keys()):
            del self.model.historyOutputRequests[request_name]

        if analysis_type == 'buckle':
            if output_mode == 'full':
                self.model.FieldOutputRequest(name='Buckling-Output',
                    createStepName=self._get_buckle_step_name(),
                    variables=buckling_variables)
            else:
                # Production ODBs retain the first buckling mode shape for
                # engineering review without storing U for every extracted
                # eigenmode.  All requested eigenvalues remain available in
                # the .dat file for the optimisation constraint.
                self.model.FieldOutputRequest(name='Buckling-Output',
                    createStepName=self._get_buckle_step_name(),
                    variables=buckling_variables, modes=(1,))
            return

        static_step_name = self._get_static_step_name()
        static_output = self.model.FieldOutputRequest(name='F-Output-1',
            createStepName=static_step_name, variables=variables,
            frequency=LAST_INCREMENT)

        if analysis_type == 'static_buckle':
            buckle_step_name = self._get_buckle_step_name()
            if output_mode == 'full':
                buckle_output = self.model.FieldOutputRequest(
                    name='Buckling-Output',
                    createStepName=buckle_step_name,
                    variables=buckling_variables)
            else:
                buckle_output = self.model.FieldOutputRequest(
                    name='Buckling-Output',
                    createStepName=buckle_step_name,
                    variables=buckling_variables, modes=(1,))

            # A field output request propagates into all subsequent steps, and
            # can only be deactivated in a step that comes after the step it
            # was created in.  The two base states use opposite step orders,
            # so deactivate whichever request is created first.
            if self._get_buckling_base_state() == 'initial':
                # Buckling -> static: do not propagate the modal request into
                # the static step; F-Output-1 supplies the displacement field
                # there.
                buckle_output.deactivate(static_step_name)
            else:
                # Static -> buckling: do not propagate the static request into
                # the perturbation step; Buckling-Output supplies the mode
                # shape there.
                static_output.deactivate(buckle_step_name)

        # Wing tip displacement for the design evaluation: history output of
        # U for the tip cross-section nodes only, a few kilobytes of data
        # that the fast post-processing reads instead of scanning the full
        # nodal U field.
        tip_set = self._create_tip_node_set()
        if tip_set is not None:
            self.model.HistoryOutputRequest(
                name='Tip-Output', createStepName=static_step_name,
                variables=('U',), region=tip_set, frequency=LAST_INCREMENT)

        # A metallic rib carries no composite layup, so it appears in none of
        # the layup requests below and its stress field is the only way to
        # check its strength.
        if rib_variables:
            a = self.rootAssembly
            for rib in self.ribs:
                if rib.name_layups:
                    continue
                self.model.FieldOutputRequest(
                    name='Rib-Output-%s' % rib.name_part,
                    createStepName=static_step_name,
                    variables=rib_variables,
                    frequency=LAST_INCREMENT,
                    region=a.instances[rib.name_part].sets['all'])

        # An empty variable list means that no ply-by-ply output is needed,
        # e.g., the optimisation mode without a failure model.
        if not layup_variables:
            print('>>> OUTPUT_VARIABLES no layup output request is created.')
            return

        layup_locations = self._get_layup_location_options()
        print('>>> OUTPUT_VARIABLES layup section points: %s'
              % self.pRun.get('layup_output_ply_locations', 'all'))

        for i, name_layup in enumerate(self.lofting.name_layups):
            self.model.FieldOutputRequest(name='Layup-Output-lofting-%d' % i,
                createStepName=static_step_name, variables=layup_variables,
                frequency=LAST_INCREMENT,
                layupNames=('lofting.%s' % name_layup, ),
                rebar=EXCLUDE, **layup_locations)

        for i, rib in enumerate(self.ribs):
            for j, name_layup in enumerate(rib.name_layups):
                self.model.FieldOutputRequest(name='Layup-Output-rib%d-%d' % (i, j),
                    createStepName=static_step_name, variables=layup_variables,
                    frequency=LAST_INCREMENT,
                    layupNames=('%s.%s' % (rib.name_part, name_layup), ),
                    rebar=EXCLUDE, **layup_locations)

    def setup_jobs(self):
        '''
        Define jobs
        '''
        mdb.Job(name=self.name_job, model=str(self.name_model), description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memory=self.pRun['memory_max_percentage'], 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, 
            numThreadsPerMpiProcess=self.pRun['numThreadsPerMpiProcess'], 
            multiprocessingMode=DEFAULT, 
            numCpus=self.pRun['numCpus'], numDomains=self.pRun['numCpus'], numGPUs=0)
        

if __name__ == '__main__':

    fname = get_parameter_file(sys.argv)
    print('>>> Reading parameters from: %s' % fname)
    parameters = load_parameters(fname)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']

    #* Build model

    name_job = str(parameters['name_job'])

    model = WingboxModel(name_job, pGeo, pMesh, pRun)
    model.build()
    model.set_view()
    model.save_cae('WingBox.cae')

    #* Total mass of the model (optimisation objective).  The composite
    #* material is created without a density card (see
    #* Model.create_material_IM785517), so assembly getMassProperties()
    #* cannot return a mass; compute it from the per-part volumes (exact
    #* shell thickness x area, including the non-design thickness factors)
    #* and the library densities instead.
    from AbaqusTools.materials import MATERIAL_LIBRARY
    density_composite = float(MATERIAL_LIBRARY[
        str(pMesh.get('material_name', 'IM7/8551-7'))]['density'])
    density_rib = float(MATERIAL_LIBRARY[
        str(pMesh.get('rib_material_name', 'Aluminum-7075'))]['density'])

    volume_of_part = {}
    mass_tonne = 0.0
    for name_part in model.model.parts.keys():
        volume = float(model.model.parts[name_part].getMassProperties()['volume'])
        density = density_rib if name_part.startswith('rib_') else density_composite
        volume_of_part[name_part] = volume
        mass_tonne += volume*density
    with open(name_job + '_mass.json', 'w') as f:
        json.dump({
            'mass_tonne': mass_tonne,
            'mass_kg': mass_tonne*1.0e3,
            'volume_mm3_of_part': volume_of_part,
            'density_composite_tonne_per_mm3': density_composite,
            'density_rib_tonne_per_mm3': density_rib,
        }, f, indent=2)
    print('>>> MASS %s: %.6f tonne (%.3f kg)' % (name_job, mass_tonne, mass_tonne*1.0e3))

    execution_mode = str(parameters.get('execution_mode', 'default')).lower()
    if execution_mode in ('build', 'build_only', 'build-only'):
        pass
    elif execution_mode in ('write_input', 'write-input'):
        model.write_job_inp()
        if get_failure_model(pMesh) == 'larc05':
            model.write_IM785517_property_table_inp(
                method=pMesh['user_subroutine'], fname_input=model.name_job+'.inp')
    elif execution_mode in ('datacheck', 'data_check', 'data-check'):
        if get_failure_model(pMesh) == 'larc05':
            raise RuntimeError(
                'In-process datacheck is not supported for the patched LaRC05 input deck.'
            )
        model.submit_job(name_job, only_data_check=True)
        job = mdb.jobs[name_job]
        running_statuses = (
            'NONE', 'SUBMITTED', 'RUNNING', 'CHECK_SUBMITTED', 'CHECK_RUNNING',
        )
        terminal_error_statuses = ('ABORTED', 'TERMINATED')
        job_log_path = name_job + '.log'
        success_marker = ('ABAQUS JOB %s COMPLETED' % name_job).upper()
        error_markers = (
            'ABAQUS/ANALYSIS EXITED WITH ERRORS',
            'ANALYSIS INPUT FILE PROCESSOR EXITED WITH AN ERROR',
        )
        wait_timeout = float(pRun.get('datacheck_timeout_seconds', 1800.0))
        wait_started = time.time()
        job_status = str(job.status).upper()
        while True:
            job_status = str(job.status).upper()
            job_log_upper = ''
            if os.path.isfile(job_log_path):
                with open(job_log_path, 'r') as job_log:
                    job_log_upper = job_log.read().upper()

            if success_marker in job_log_upper:
                job_status = 'CHECK_COMPLETED'
                break
            if any(marker in job_log_upper for marker in error_markers):
                if job_status in running_statuses:
                    job_status = 'ABORTED'
                break
            if job_status == 'CHECK_COMPLETED':
                break
            if job_status in terminal_error_statuses:
                break
            if time.time() - wait_started > wait_timeout:
                raise RuntimeError(
                    'Timed out after %.1f seconds waiting for datacheck %s. '
                    'Last CAE status: %s; job log: %s.'
                    % (wait_timeout, name_job, job_status, job_log_path)
                )
            time.sleep(1.0)
        print('>>> DATACHECK_STATUS %s %s' % (name_job, job_status))
        sys.stdout.flush()
        if job_status != 'CHECK_COMPLETED':
            raise RuntimeError(
                'Abaqus datacheck failed for %s with status %s.'
                % (name_job, job_status)
            )
    elif execution_mode not in ('default', 'auto'):
        raise ValueError('Invalid execution_mode: %s' % execution_mode)
    elif get_failure_model(pMesh) == 'larc05':
        model.write_job_inp()
        model.write_IM785517_property_table_inp(
            method=pMesh['user_subroutine'], fname_input=model.name_job+'.inp')
    else:
        if not parameters['not_run_job']:
            model.submit_job(name_job)
            # Abaqus/CAE submits asynchronously.  Returning here would end the
            # CAE session while the solver is still starting up, and the caller
            # would clean away the *.com/*.env/*.sim/*.stt files it needs, so
            # waiting is the default.
            if bool(parameters.get('wait_for_completion', True)):
                print('>>> Waiting for Abaqus job %s to complete.' % name_job)
                sys.stdout.flush()
                mdb.jobs[name_job].waitForCompletion()
                print('>>> ABAQUS_SOLVE_STATUS %s %s'
                      % (name_job, str(mdb.jobs[name_job].status)))
                sys.stdout.flush()
    
