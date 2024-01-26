# class ParallelPack:
    # def __init__(self, num_cells: int, model: pybamm.BaseModel, geo: dict):
        # self.model = model
        # self.geo = geo

        # self.num_cells = num_cells
        # self.cells = [Cell(f"Cell{i + 1}", model, geo) for i in range(NUM_CELLS)]

        # self.iapp_total_param = pybamm.Parameter("Input Current / Area")

        # # Vcell1 - Vcell2 = 0
        # # (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
        # model.algebraic = {}
        # for i in range(len(self.cells) - 1):
            # model.algebraic[self.cells[i].voltage] = self.cells[i + 1].voltage - self.cells[i].voltage
            # model.algebraic[self.cells[i].iapp] = self.cells[i].pos.phi - self.cells[i].neg.phi - self.cells[i].voltage

        # model.algebraic[self.cells[-1].voltage] = self.cells[-1].pos.phi - self.cells[-1].neg.phi - self.cells[-1].voltage
        # model.algebraic[self.cells[-1].iapp] = i_total - sum(cell.iapp for cell in self.cells)