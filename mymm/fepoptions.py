import logging

class FepOptions:
    def __init__(self, output_file,
                 num_initial_equil_steps,
                 num_initial_min_steps,
                 num_total_steps_per_window,
                 num_equil_steps_per_window,
                 prefix = "",
                 output_freq = 10,
                 alch_vdw_lambda_end = 0.0,
                 alch_elec_lambda_start = 0.0,
                 alch_vdw_shift_coeff = 5.0,
                 alch_decouple = "on"):
        self.prefix                 = prefix
        self.num_initial_equil_steps= num_initial_equil_steps
        self.num_initial_min_steps  = num_initial_min_steps
        self.output_file            = output_file
        self.output_freq            = output_freq
        self.alch_vdw_lambda_end    = alch_vdw_lambda_end
        self.alch_elec_lambda_start = alch_elec_lambda_start
        self.alch_vdw_shift_coeff   = alch_vdw_shift_coeff
        self.alch_decouple          = alch_decouple

        if num_total_steps_per_window < num_equil_steps_per_window:
            print("Error! num_total_steps_per_window < num_equil_steps_per_window!\n")
            return
        self.alch_equil_steps           = num_equil_steps_per_window
        self.num_total_steps_per_window = num_total_steps_per_window

    def __str__(self):
        mystr = ""

        mystr += "\tprefix = %s\n"%(self.prefix)
        mystr += "\toutput_file = %s\n"%(self.output_file)
        mystr += "\tnum_initial_equil_steps = %d\n"%(self.num_initial_equil_steps)
        mystr += "\tnum_initial_min_steps   = %d\n"%(self.num_initial_min_steps)
        mystr += "\tnum_equil_steps_per_window = %d\n"%(self.alch_equil_steps)
        mystr += "\tnum_total_steps_per_window = %d\n"%(self.num_total_steps_per_window)
        mystr += "\toutput_freq = %d\n"%(self.output_freq)
        mystr += "\talch_vdw_lambda_end = %f\n"%(self.alch_vdw_lambda_end)
        mystr += "\talch_elec_lambda_start = %f\n"%(self.alch_elec_lambda_start)
        mystr += "\talch_vdw_shift_coeff = %f\n"%(self.alch_vdw_shift_coeff)
        mystr += "\talch_decouple = %s\n"%(self.alch_decouple)
        return mystr

    def get_equilibration_commands(self, temp, lambda_point = 0.0):
        return ["runFEPmin %f %f 0.0 %d %d %f"%(lambda_point, lambda_point, self.num_initial_equil_steps, self.num_initial_min_steps, temp)]

    def get_production_commands(self, template_command_list):
        command_list = []

        for command in template_command_list:
            command_list.append(command + "  %d"%(self.num_total_steps_per_window))

        return command_list
