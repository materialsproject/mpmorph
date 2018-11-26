from fireworks import explicit_serialize, FireTaskBase, FWAction, Workflow
from pymatgen import Structure


@explicit_serialize
class DiffusionTask(FireTaskBase):
    required_params = ['temperatures', 'max_steps', 'target_steps', 'trajectory_to_db']
    optional_params = []

    def run_task(self, fw_spec):
        s = Structure.from_file('CONTCAR.gz')
        fws = []
        for t in self['temperatures']:
            fws.extend(get_converge_new(s, t, max_steps=self['max_steps'],
                                           target_steps=self['target_steps'],
                                           trajectory_to_db=self['trajectory_to_db']))
        wf = Workflow(fws)
        return FWAction(detours=wf)


from mpmorph.workflow.converge import get_converge_new