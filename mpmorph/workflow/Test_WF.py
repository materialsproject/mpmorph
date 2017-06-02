from fireworks import explicit_serialize, FireTaskBase, FWAction, Firework, LaunchPad, Workflow

def get_dummy_wf():
    priority_spec = {"_priority": 100000}
    t = EmptyTask()

    fws = []
    fws.append(Firework(tasks=[t], name="run0", parents=[], spec=priority_spec))

    t = SpawnTask(spawn_count=0, spawn_limit=3, priority_spec=priority_spec)
    fws.append(Firework(tasks=[t], name="spawn_0", parents=fws[len(fws)-1], spec=priority_spec))

    wf = Workflow(fireworks=fws, name="Empty_WF_01")
    return wf

@explicit_serialize
class EmptyTask(FireTaskBase):

    def run_task(self, fw_spec):
        return FWAction()

@explicit_serialize
class SpawnTask(FireTaskBase):
    required_params = ["spawn_count", "spawn_limit", "priority_spec"]
    optional_params = ["diffusion"]
    def run_task(self, fw_spec):
        spawn_count = self["spawn_count"]
        spawn_limit = self["spawn_limit"]
        priority_spec = self["priority_spec"]
        diffusion = self.get("diffusion", False)
        fws = []
        if spawn_count < spawn_limit:
            fws.append(Firework(SpawnTask(spawn_count=spawn_count+1, spawn_limit=spawn_limit, priority_spec=priority_spec, diffusion=diffusion),
                                name="spawn_"+ str(spawn_count+1), spec=priority_spec))
        else:
            if diffusion:
                fws.append(Firework(EmptyTask(), name="longrun_0", spec=priority_spec))
                fws.append(Firework(EmptyTask(), name="longrun_1", spec=priority_spec, parents=fws[len(fws)-1]))
                fws.append(Firework(EmptyTask(), name="longrun_2", spec=priority_spec, parents=fws[len(fws)-1]))
                fws.append(Firework(EmptyTask(), name="longrun_3", spec=priority_spec, parents=fws[len(fws)-1]))
                fws.append(Firework(EmptyTask(), name="longrun_4", spec=priority_spec, parents=fws[len(fws)-1]))
            else:
                fws.append(Firework(EmptyTask(), name="longrun_0", spec=priority_spec))
                fws.append(Firework(EmptyTask(), name="longrun_1", spec=priority_spec, parents=fws[len(fws)-1]))
                t=SamplerTask(spawn_count=10, priority_spec=priority_spec)
                fws.append(Firework(t, name="sampler", spec=priority_spec, parents=fws[len(fws)-1]))
        wf = Workflow(fws)
        return FWAction(additions=wf)

@explicit_serialize
class SamplerTask(FireTaskBase):
    required_params = ["spawn_count", "priority_spec"]
    optional_params = []
    def run_task(self, fw_spec):
        spawn_count=self["spawn_count"]
        priority_spec = self["priority_spec"]
        fws=[]
        for i in range(spawn_count):
            fws.append(Firework(EmptyTask(), name="snap_" + str(i) + "_cool_1500", spec=priority_spec))
            fws.append(Firework(EmptyTask(), name="snap_" + str(i) + "_hold_1500", spec=priority_spec, parents=fws[len(fws)-1]))
            fws.append(Firework(EmptyTask(), name="snap_" + str(i) + "_cool_1000", spec=priority_spec, parents=fws[len(fws)-1]))
            fws.append(Firework(EmptyTask(), name="snap_" + str(i) + "_hold_1000", spec=priority_spec, parents=fws[len(fws)-1]))
            final_fw = len(fws)-1
            fws.append(Firework(EmptyTask(), name="relax_0", spec=priority_spec, parents=fws[final_fw]))
            fws.append(Firework(EmptyTask(), name="relax_1", spec=priority_spec, parents=fws[len(fws)-1]))
            fws.append(Firework(EmptyTask(), name="relax_2", spec=priority_spec, parents=fws[len(fws)-1]))
            if i==0:
                for temp in [500, 1000, 1500]:
                    fws.append(Firework(EmptyTask(), name="diffusion_run0_" + str(temp), spec=priority_spec, parents=fws[final_fw]))
                    fws.append(Firework(SpawnTask(spawn_count=0, spawn_limit=3, diffusion=True, priority_spec=priority_spec),
                                        name="spawn_1", spec=priority_spec, parents=fws[len(fws)-1]))
        wf = Workflow(fws)
        return FWAction(additions=wf)