#%%
import os
from pathlib import Path
import pandas as pd
from bids2table import bids2table

#%%
class getfiles:
    def __init__(self, BIDS_path):
        self.BIDS_path = BIDS_path
        self.BIDS_df = bids2table(self.BIDS_path)

    def locT1w(self):
        """
        Get the location of T1w images for all subjects in the BIDS dataset
        """
        bidsdf = self.BIDS_df.filter_multi(
            suffix="T1w",
            ses={"items": ["implant01", "preimplant01"]}
        )
        bidsdf = bidsdf.flat

        # drop the columns that are all missing
        bidsdf = bidsdf.dropna(axis=1, how='all')

        t1w = bidsdf.filter(items=['sub', 'ses', 'file_path'])
        t1w.set_index('sub', inplace=True)

        return t1w
    
    def locCT(self):
        """
        Get the location of CT images for all subjects in the BIDS dataset
        """
        bidsdf = self.BIDS_df.filter_multi(
            suffix={"items": ["CT", "ct"]},
            ses={"items": ["implant01", "preimplant01"]}
        )
        bidsdf = bidsdf.flat

        # drop the columns that are all missing
        bidsdf = bidsdf.dropna(axis=1, how='all')

        ct = bidsdf.filter(items=['sub', 'ses', 'file_path'])
        ct.set_index('sub', inplace=True)

        return ct
    
    def electrodes(self):
        """
        Get the location of electrodes for all subjects in the BIDS dataset
        """
        bidsdf = self.BIDS_df.filter_multi(
            suffix="electrodes",
            ses={"items": ["implant01", "preimplant01"]}
        )
        bidsdf = bidsdf.flat

        # drop the columns that are all missing
        bidsdf = bidsdf.dropna(axis=1, how='all')

        electrodes = bidsdf.filter(items=['sub', 'ses', 'file_path'])
        electrodes.set_index('sub', inplace=True)
        

        return electrodes

#%%
if __name__ == "__main__":
    project_path = Path(os.getenv('PROJECT_ROOT', os.getcwd()))
    BIDS_path = project_path / 'data'
    files = getfiles(BIDS_path)
    t1 = files.locT1w()
    ct = files.locCT()
    electrodes = files.electrodes()

    file_paths = pd.concat([
        t1['file_path'].rename('t1w'),
        ct['file_path'].rename('ct'),
        electrodes['file_path'].rename('electrodes')
    ], axis=1).to_csv(project_path / 'code' / 'file_paths.csv')
#%%
