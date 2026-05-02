import ROOT as r

RNTupleModel = r.RNTupleModel
RNTupleReader = r.RNTupleReader
RNTupleWriter = r.RNTupleWriter

class RNtupleMaker:
    def __init__(self, ofile, columns_add):
        self.fname = ofile.split(":")[0]
        self.tname = ofile.split(":")[1]
        self.model = RNTupleModel.Create()
        for vartype, list_names in columns_add.items():
            for name in list_names:
                print(f"adding type = {vartype} with name = {name}")
                self.model.MakeField[vartype](f"{name}")
        self.writer = RNTupleWriter.Recreate(self.model, self.tname,
                                             self.fname)
