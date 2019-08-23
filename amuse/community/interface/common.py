"""
Common code functions
"""

from amuse.support.interface import InCodeComponentImplementation

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification

import sys
class CommonCodeInterface(object):

    @legacy_function
    def initialize_code():
        """
        Run the initialization for the code, called before
        any other call on the code (so before any parameters
        are set or particles are defined in the code).
        """

        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
        -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        -2 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def cleanup_code():
        """
        Run the cleanup for the code, called
        just before stopping the code. No functions
        should be called after this code.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
        -1 - ERROR
            Error happened during cleanup, this error needs to be further specified by every code implemention
        -2 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def commit_parameters():
        """
        Perform initialization in the code dependent on the
        values of the parameters.
        Called after the parameters have been set or updated.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
        -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def recommit_parameters():
        """
        Perform initialization actions after parameters
        have been updated (after commit_parameters and
        particles have been loaded).
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """

        return function

    def invoke_state_change(self):
        pass
        
        

class CommonCode(InCodeComponentImplementation):

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('INITIALIZED', 'before_get_parameter')
        object.add_method('INITIALIZED', 'before_set_parameter')
        object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')
    
    def define_methods(self, object):
        object.add_method(
            'initialize_code',
            (),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'cleanup_code',
            (),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'commit_parameters',
            (),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'recommit_parameters',
            (),
            (object.ERROR_CODE)
        )
    
