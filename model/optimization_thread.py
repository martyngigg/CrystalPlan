import threading
from model.optimization import run_optimization
from model.utils import cleanup_wxpython

class OptimizationThread(threading.Thread):
    """Thread to run the GA optimization."""
    def __init__ (self, controller):
        threading.Thread.__init__(self)
        # self.controller = controller
        self.controller = controller
        self.start() #Start on creation

    def run(self):
        #Just run the optimization
        self.controller.params.optimization_running = True
        try:
            (ga, aborted, converged) = run_optimization(self.controller.params, self.controller.step_callback)
            # self.controller.params.optimization_running = False
            #Call the completion function.
            self.controller.complete( ga, aborted, converged )
        except Exception as e:
            print "Error while running optimization", e
            self.controller.restore_buttons()
        finally:
            pass
            self.controller.params.optimization_running = False

        cleanup_wxpython()