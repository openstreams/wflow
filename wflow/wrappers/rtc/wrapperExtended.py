import bmi
import bmi.wrapper


class BMIWrapperExtended(bmi.wrapper.BMIWrapper):

    # base
    def update_until(self, time):

        t = self.get_current_time()
        tEnd = self.get_end_time()
        while t < time:
            self.update()
            t = self.get_current_time()

    # getter / setter
    def get_value(self, var_name):
        return self.get_var(var_name)

    def set_value(self, var_name, src):
        self.set_var(var_name, src)

    # info
    def get_component_name(self):
        return "RTC-Tools"

    # time
    def get_time_units(self):
        return "s"

    # vars
    def get_var_itemsize(self, var_name):
        return 8

    def get_var_nbytes(self, var_name):
        return 8
