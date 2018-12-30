from ceci import PipelineStage
from .types import FitsFile

class ReduceCat(PipelineStage) :
    name="ReduceCat"
    inputs=[('raw_data',FitsFile)]
    outputs=[('clean_catalog',FitsFile)]
    #outputs=[('clean_catalog',FitsFile),('masked_fraction',FitsFile),
    #         ('depth_map',FitsFile),('star_map',FitsFile)]
    config_options={'min_snr':10.,'sv_depth':True,'sv_mask':True,
                    'sv_syst':True,'depth_cut':24.5,'res':0.0285,
                    'res_bo':0.003,'pad':0.1,'band':'i','depth_method':2,
                    'flat_project':'CAR','mask_type':'sirius'}

    def run(self) :
        for inp,_ in self.inputs :
            fname=self.get_input(inp)
            print("Reading "+fname)
            open(fname)

        for out,_ in self.outputs :
            fname=self.get_output(out)
            print("Writing "+fname)
            open(fname,"w")        
                    
if __name__ == '__main__':
    cls = PipelineStage.main()
