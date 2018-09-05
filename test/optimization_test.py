import unittest

import numpy as np
import random
import model.experiment as experiment
import model.instrument as instrument
import model.goniometer as goniometer
from model.optimization import run_optimization, OptimizationParameters
from model.instrument import PositionCoverage

class OptimizationTests(unittest.TestCase):

    def get_default_parameters(self):
        op = OptimizationParameters()
        op.desired_coverage = 85
        op.number_of_orientations = 10
        op.max_generations = 10
        op.population = 10
        op.use_multiprocessing = False
        op.use_volume = False
        return op
    
    def setUp(self):
        np.random.seed(42)
        random.seed(42)

        instrument.inst = instrument.Instrument("instruments/TOPAZ_geom_all_2011.csv")
        instrument.inst.set_goniometer(goniometer.TopazAmbientGoniometer())
        instrument.inst.positions = [PositionCoverage([0.0, 0.0], None, np.identity(3)), PositionCoverage([1.0, 0.0], None, np.identity(3))]
                                                  
        experiment.exp = experiment.Experiment(instrument.inst)
        experiment.exp.initialize_reflections()
        print('setup')
        experiment.exp.recalculate_reflections(None)
        experiment.exp.verbose = False

    def test_run_optimization_single_process(self):
        op = self.get_default_parameters()

        (ga, aborted, converged) = run_optimization(op)

        self.assertFalse(aborted)
        self.assertFalse(converged)

        best = ga.bestIndividual()
        self.assertAlmostEqual(best.coverage*100.0, 43.197, places=2)
        assert False
        
    def test_run_optimization_multi_process(self):
        op = self.get_default_parameters()
        op.use_multiprocessing = True
        op.number_of_processors = 4

        (ga, aborted, converged) = run_optimization(op)

        self.assertFalse(aborted)
        self.assertFalse(converged)

        best = ga.bestIndividual()
        self.assertAlmostEqual(best.coverage*100.0, 43.197, places=2)
        assert False