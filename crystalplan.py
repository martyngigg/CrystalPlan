#!/usr/bin/env python
"""The CrystalPlan application.

CrystalPlan is an experiment planning tool for crystallography.
You can choose an instrument and supply your sample's lattice
parameters to simulate which reflections will be measured,
by which detectors and at what wavelengths.

Author: Janik Zikovsky, zikovskyjl@ornl.gov
Version: $Id$
"""
# Author: Janik Zikovsky, zikovskyjl@ornl.gov
# Version: $Id$

#Simply import and launch the GUI
from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'wx'

import vtk
errOut = vtk.vtkFileOutputWindow()
errOut.SetFileName("VTK Error Out.txt")
vtkStdErrOut = vtk.vtkOutputWindow()
vtkStdErrOut.SetInstance(errOut)

import gui.main
import multiprocessing
multiprocessing.freeze_support()

gui.main.handle_arguments_and_launch(InstalledVersion=True)
