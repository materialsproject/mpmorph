from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction

@explicit_serialize
class ConvergeTask(FireTaskBase):
    """

    Ensures a structure is converged before production MD run

    """

    required_params = ["conv_parameter"]
    optional_params = []

    def run_task(self, fw_spec):
        vasp_cmd = self["vasp_cmd"]