from ...core.base_components import Monitor
import matplotlib.pyplot as plt


class PlotFunctionMonitor(Monitor):

    def __init__(self, plot_function):
        plt.ion()
        self._plot_function = plot_function
        self._fig = plt.figure()

    def store(self, state):
        """
        Stores the given state in the Monitor and performs class-specific
        actions.

        Args:
            state (dict): A model state dictionary.
        """
        self._fig.clear()
        self._plot_function(self._fig, state)
        plt.draw_all()
        plt.pause(1e-5)  # necessary to draw, pause can be arbitrarily small
