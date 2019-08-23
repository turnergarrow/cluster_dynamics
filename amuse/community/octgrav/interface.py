from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class OctgravInterface(
    CodeInterface, 
    LiteratureReferencesMixIn,
    GravitationalDynamicsInterface,
    StoppingConditionInterface,
    GravityFieldInterface):
    """
        .. [#] Gaburov E., Bedorf J., Portegies Zwart S., "A gravitational tree code on graphics processing units:
               Implementation in CUDA", 2010, Proc. C. Sc., 1, 1119; and main MUSE paper, arXiv/0807.1996
    """

    include_headers = ['interface.h', 'parameters.h', 'worker_code.h', 'local.h', 'stopcond.h']

    def __init__(self, convert_nbody = None, **options):
        CodeInterface.__init__(self, name_of_the_worker="octgrav_worker", **options)
        """
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()

        self.convert_nbody = convert_nbody
        """
        LiteratureReferencesMixIn.__init__(self)

    

    

    @legacy_function
    def get_theta_for_tree():
        """
        
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('opening_angle', dtype='float64', direction=function.OUT,
            description = "opening_angle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            xx
        -1 - ERROR
            xx
        """
        return function    

    @legacy_function
    def set_theta_for_tree():
        """
        
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('opening_angle', dtype='float64', direction=function.IN,
            description = "opening_angle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            xx
        -1 - ERROR
            xx
        """
        return function    

    @legacy_function
    def set_time_step():
        """
        Update timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function
        


        


class Octgrav(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody = None, **options):
        legacy_interface = OctgravInterface(**options)
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations", 
            default_value = 0.01 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = 0.01 | nbody_system.time
        )
        object.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle",
            "opening angle for building the tree between 0 and 1", 
            default_value = 0.8
        )
        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
        
        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_step",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_theta_for_tree",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_theta_for_tree",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(object)


    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)
        
# this should be checked!
        object.add_method('EDIT', 'get_gravity_at_point')
        object.add_method('EDIT', 'get_potential_at_point')
        
        self.stopping_conditions.define_state(object)
        
        
    

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)
