import unittest

import numpy as np
import random
import model.experiment as experiment
import model.instrument as instrument
import model.goniometer as goniometer
from model.optimization import run_optimization, OptimizationParameters
from model.optimization_thread import OptimizationThread
from model.instrument import PositionCoverage
from gui.frame_optimizer import OptimizerController
import wx

# Helper functions for tests
# ------------------------------------------------------------------------------------------
def get_default_parameters():
    op = OptimizationParameters()
    op.desired_coverage = 85
    op.number_of_orientations = 5
    op.max_generations = 2
    op.population = 20
    op.use_multiprocessing = False
    op.use_volume = False
    return op

def setup_instrument():
    instrument.inst = instrument.Instrument("instruments/TOPAZ_geom_all_2011.csv")
    instrument.inst.set_goniometer(goniometer.TopazAmbientGoniometer())
    instrument.inst.verbose = True

def setup_experiment():
    experiment.exp = experiment.Experiment(instrument.inst)
    experiment.exp.initialize_reflections()
    experiment.exp.recalculate_reflections(None)
    experiment.exp.verbose = False

# ------------------------------------------------------------------------------------------

class MockOptimizerController:
    completed_callback_called = False

    def __init__(self, params):
        self.params = params
        self._want_abort = False
    
    def step_callback(self, ga):
        pass

    def complete(self, ga, abort, converged):
        self.abort = abort
        self.converged = converged
        self.completed_callback_called = True
    
    def restore_buttons(self):
        pass

class OptimizationThreadTests(unittest.TestCase):

    def setUp(self):
        np.random.seed(42)
        random.seed(42)

        setup_instrument()
        setup_experiment()

    def test_run_optimization_thread(self):
        app = wx.App()

        params = get_default_parameters()
        params.use_multiprocessing = True
        params.number_of_processors = 4
        
        controller = MockOptimizerController(params)

        thread = OptimizationThread(controller)
        thread.join()

        self.assertFalse(controller.abort)
        self.assertFalse(controller.converged)
        self.assertTrue(controller.completed_callback_called)

# ------------------------------------------------------------------------------------------

class OptimizationTests(unittest.TestCase):
    
    def setUp(self):
        np.random.seed(42)
        random.seed(42)

        setup_instrument()
        setup_experiment()

    def test_run_optimization_single_process(self):
        op = get_default_parameters()

        (ga, aborted, converged) = run_optimization(op)

        self.assertFalse(aborted)
        self.assertFalse(converged)

        best = ga.bestIndividual()
        self.assertAlmostEqual(best.coverage*100.0, 27.079, places=2)
    
    def test_run_optimization_multi_process(self):
        op = get_default_parameters()
        op.use_multiprocessing = True
        op.number_of_processors = 4

        (ga, aborted, converged) = run_optimization(op)

        self.assertFalse(aborted)
        self.assertFalse(converged)

        best = ga.bestIndividual()
        self.assertAlmostEqual(best.coverage*100.0, 27.079, places=2)

if __name__ == "__main__":
    unittest.main()