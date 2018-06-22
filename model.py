import modelbase


class Poolman(modelbase.Model):
    def __init__(self, p, r):
        super().__init__(p)
        compounds = ["PGA", "BPGA", "GAP", "DHAP", "FBP", "F6P", "G6P", "G1P", "SBP", "S7P",
                     "E4P", "X5P", "R5P", "RUBP", "RU5P", "ATP"]
        self.set_cpds(compounds)
        self.add_algebraicModule(r.ADP, "ADP_mod", ["ATP"], ["ADP"])
        self.add_algebraicModule(r.P_i, "Pi_mod", ["PGA", "BPGA", "GAP", "DHAP", "FBP", "F6P", "G6P", "G1P", "SBP", "S7P", "E4P", "X5P", "R5P", "RUBP", "RU5P", "ATP"], ["P"])
        self.add_algebraicModule(r.N, "N_mod", ["P", "PGA", "GAP", "DHAP"], ["N"])
        self.add_reaction("v1", r.v1, {"PGA": 2, "RUBP": -1}, "RUBP", "PGA", "FBP", "SBP", "P")
        self.add_reaction("v2", r.v2, {"PGA": -1, "ATP": -1, "BPGA": 1}, "ATP", "PGA", "ADP", "BPGA")
        self.add_reaction("v3", r.v3, {"BPGA": -1, "GAP": 1}, "BPGA", "GAP", "P")
        self.add_reaction("v4", r.v4, {"GAP": -1, "DHAP": 1}, "GAP", "DHAP")
        self.add_reaction("v5", r.v5, {"GAP": -1, "DHAP": -1, "FBP": 1}, "GAP", "DHAP", "FBP")
        self.add_reaction("v6", r.v6, {"FBP": -1, "F6P": 1}, "FBP", "F6P", "P")
        self.add_reaction("v7", r.v7, {"GAP": -1, "F6P": -1, "E4P": 1, "X5P": 1}, "GAP", "F6P", "X5P", "E4P")
        self.add_reaction("v8", r.v8, {"DHAP": -1, "E4P": -1, "SBP": 1}, "DHAP", "E4P", "SBP")
        self.add_reaction("v9", r.v9, {"SBP": -1, "S7P": 1}, "SBP", "P")
        self.add_reaction("v10", r.v10, {"GAP": -1, "S7P": -1, "X5P": 1, "R5P": 1}, "GAP", "S7P", "X5P", "R5P")
        self.add_reaction("v11", r.v11, {"R5P": -1, "RU5P": 1}, "R5P", "RU5P")
        self.add_reaction("v12", r.v12, {"X5P": -1, "RU5P": 1}, "X5P", "RU5P")
        self.add_reaction("v13", r.v13, {"ATP": -1, "RU5P": -1, "RUBP": 1}, "RU5P", "ATP", "PGA", "RUBP", "P", "ADP")
        self.add_reaction("v14", r.v14, {"F6P": -1, "G6P": 1}, "F6P", "G6P")
        self.add_reaction("v15", r.v15, {"G6P": -1, "G1P": 1}, "G6P", "G1P")
        self.add_reaction("v16", r.v16, {"ATP": 1}, "ADP", "P")
        self.add_reaction("vPGA_out", r.vPGA_out, {"PGA": -1}, "PGA", "N")
        self.add_reaction("vGAP_out", r.vGAP_out, {"GAP": -1}, "GAP", "N")
        self.add_reaction("vDHAP_out", r.vDHAP_out, {"DHAP": -1}, "DHAP", "N")
        self.add_reaction("vSt", r.vStarchProduction, {"G1P": -1, "ATP": -1},
                          "G1P", "ATP", "ADP", "P", "PGA", "F6P", "FBP")


class PoolmanNADPH(modelbase.Model):
    def __init__(self, p, r):
        super().__init__(p)
        compounds = ["PGA", "BPGA", "GAP", "DHAP", "FBP", "F6P", "G6P", "G1P", "SBP", "S7P",
                     "E4P", "X5P", "R5P", "RUBP", "RU5P", "ATP", "NADPH"]
        self.set_cpds(compounds)
        self.add_algebraicModule(r.ADP, "ADP_mod", ["ATP"], ["ADP"])
        self.add_algebraicModule(r.NADP, "NADPH_mod", ["NADPH"], ["NADP"])
        self.add_algebraicModule(r.P_i, "Pi_mod", ["PGA", "BPGA", "GAP", "DHAP", "FBP", "F6P", "G6P", "G1P", "SBP",
                                                   "S7P", "E4P", "X5P", "R5P", "RUBP", "RU5P", "ATP"], ["P"])
        self.add_algebraicModule(r.N, "N_mod", ["P", "PGA", "GAP", "DHAP"], ["N"])
        self.add_reaction("v1", r.v1, {"PGA": 2, "RUBP": -1}, "RUBP", "PGA", "FBP", "SBP", "P", "NADPH")
        self.add_reaction("v2", r.v2, {"PGA": -1, "ATP": -1, "BPGA": 1}, "ATP", "PGA", "ADP", "BPGA")
        self.add_reaction("v3", r.v3, {"BPGA": -1, "NADPH": -1, "GAP": 1}, "NADPH", "BPGA", "GAP", "NADP", "P")
        self.add_reaction("v4", r.v4, {"GAP": -1, "DHAP": 1}, "GAP", "DHAP")
        self.add_reaction("v5", r.v5, {"GAP": -1, "DHAP": -1, "FBP": 1}, "GAP", "DHAP", "FBP")
        self.add_reaction("v6", r.v6, {"FBP": -1, "F6P": 1}, "FBP", "F6P", "P")
        self.add_reaction("v7", r.v7, {"GAP": -1, "F6P": -1, "E4P": 1, "X5P": 1}, "GAP", "F6P", "X5P", "E4P")
        self.add_reaction("v8", r.v8, {"DHAP": -1, "E4P": -1, "SBP": 1}, "DHAP", "E4P", "SBP")
        self.add_reaction("v9", r.v9, {"SBP": -1, "S7P": 1}, "SBP", "P")
        self.add_reaction("v10", r.v10, {"GAP": -1, "S7P": -1, "X5P": 1, "R5P": 1}, "GAP", "S7P", "X5P", "R5P")
        self.add_reaction("v11", r.v11, {"R5P": -1, "RU5P": 1}, "R5P", "RU5P")
        self.add_reaction("v12", r.v12, {"X5P": -1, "RU5P": 1}, "X5P", "RU5P")
        self.add_reaction("v13", r.v13, {"ATP": -1, "RU5P": -1, "RUBP": 1}, "RU5P", "ATP", "PGA", "RUBP", "P", "ADP")
        self.add_reaction("v14", r.v14, {"F6P": -1, "G6P": 1}, "F6P", "G6P")
        self.add_reaction("v15", r.v15, {"G6P": -1, "G1P": 1}, "G6P", "G1P")
        self.add_reaction("v16", r.v16, {"ATP": 1}, "ADP", "P")
        self.add_reaction("vPGA_out", r.vPGA_out, {"PGA": -1}, "PGA", "N")
        self.add_reaction("vGAP_out", r.vGAP_out, {"GAP": -1}, "GAP", "N")
        self.add_reaction("vDHAP_out", r.vDHAP_out, {"DHAP": -1}, "DHAP", "N")
        self.add_reaction("vSt", r.vStarchProduction, {"G1P": -1, "ATP": -1},
                          "G1P", "ATP", "ADP", "P", "PGA", "F6P", "FBP")
        self.add_reaction("vNADPH", r.vNADPH, {"NADPH": 1}, "NADP")
