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
from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'wx'
#Simply import and launch the GUI
import gui.main
gui.main.handle_arguments_and_launch(InstalledVersion=True)
