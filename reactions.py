import numpy as np

class Reactions:
    ###########################################################################
    # Algebraic
    ###########################################################################
    def ADP(self, p, ATP):
        return p.Ca-ATP

    def P_i(self, p, y):
        PGA, BPGA, GAP, DHAP, FBP, F6P, G6P, G1P, SBP, S7P, E4P, X5P, R5P, RUBP, RU5P, ATP = y
        return p.Cp - (PGA + 2*BPGA + GAP + DHAP + 2*FBP + F6P + G6P + G1P
                       + 2*SBP + S7P + E4P + X5P + R5P + 2*RUBP + RU5P + ATP)

    def N(self, p, y):
        """Used to calculate vPGA_out, vGAP_out and vDHAP_out"""
        P, PGA, GAP, DHAP = y
        return (1+(1+(p.Kpxt/p.Pext))*((P/p.Kpi)+(PGA/p.Kpga)+(GAP/p.Kgap)+(DHAP/p.Kdhap)))

    ###########################################################################
    # Reaction rates
    ###########################################################################
    def v1(self, p, RUBP, PGA, FBP, SBP, P):
        """ Irreversible reaction: RuBP + CO2 -> PGA

        3 Ribulose-1, 5-bisphosphate + 3 CO2
        -- RuBisCO -->
        6 3-Phosphoglycerate
        """
        return (p.V1*RUBP)/(RUBP+p.Km1*(1+(PGA/p.Ki11)+(FBP/p.Ki12)+(SBP/p.Ki13)+(P/p.Ki14)+(p.NADPH/p.Ki15)))

    def v2(self, p, ATP, PGA, ADP, BPGA):
        """ Reversible, fast EQ reaction : PGA + ATP <-> BPGD + ADP

        6 3-Phosphoglycerate + 6 ATP
        -- Phosphoglycerate kinase (PGK) -->
        6 1, 3-Bisphosphoglycerate + 6 ADP
        """
        return p.kRE*((ATP*PGA) - (ADP*BPGA)/p.q2)

    def v3(self, p, BPGA, GAP, Phosphate_i):
        """ Reversible, fast EQ reaction: BPGA + NADPH <-> GAP + NADP

        6 1, 3-Bisphosphoglycerate + 6 NADPH + 6 H+
        -- Glyceraldehyde 3-phosphate dehydrogenase (GADPH)-->
        1 G3P + 5 Glyceraldehyde 3-phosphate

        Stroma pH is assumed to be constant
        """
        return p.kRE*((p.NADPH*BPGA*p.protonsStroma) - (1/p.q3)*(GAP*p.NADP*Phosphate_i))

    def v4(self, p, GAP, DHAP):
        """ Reversible, fast EQ reaction: GAP <-> DHAP

        (5) Glyceraldehyde 3-phosphate
        -- Triose phosphate isomerae (TPI)-->
        (?) Dihydroxyacetone phosphate
        """
        return p.kRE*((GAP) - (DHAP)/p.q4)

    def v5(self, p, GAP, DHAP, FBP):
        """ Reversible, fast EQ reaction: GAP + DHAP <-> FBP

        (5) Glyceraldehyde 3-phosphate + (?) Dihydroxyacetone phosphate
        -- Aldolase (ALD)-->
        Frucose 1, 6-bisphosphate
        """
        return p.kRE*((GAP*DHAP) - (FBP)/p.q5)

    def v6(self, p, FBP, F6P, P):
        """ Irreversible reaction: FBP -> F6P

        (?) Fructose 1, 6-bisphosphate + (?) H20
        --Fructose 1, 6-bisphosphatase (FBPase) -->
        (6) Fructose 6-phosphate + (Pi)
        """
        return (p.V6*FBP)/(FBP+p.Km6*(1+(F6P/p.Ki61)+(P/p.Ki62)))

    def v7(self, p, GAP, F6P, X5P, E4P):
        """ Reversible, fast EQ reaction: GAP + F6P <-> X5P + E4P

        (?) Fructose 6-phosphate + (?) Glyceraldehyde 3-phosphate
        -- Transketolase (TK) -->
        (?) Xylulose 5-phosphate + (?) Erythrose 4-phosphate
        """
        return p.kRE*((GAP*F6P) - (X5P*E4P)/p.q7)

    def v8(self, p, DHAP, E4P, SBP):
        """ Reversible, fast EQ reaction: DHAP + E4P <-> SBP

        (?) Dihydroxyacetone phosphate + (?) Erythrose 4-phosphate
        -- Aldolase (ALD)-->
        (?) Sedoheptulose 1, 7-bisphosphate
        """
        return p.kRE*((DHAP*E4P) - (SBP)/p.q8)

    def v9(self, p, SBP, P):
        """ Irreversible reaction: SBP -> S7P

        (?) Sedoheptulose 1, 7-bisphosphate + H20
        --Sedoheptulose 1, 7-bisphosphatase (SBPase)-->
        (?) Sedoheptulose 7-phosphate + (?) Pi
        """
        return (p.V9*SBP)/(SBP+p.Km9*(1+(P/p.Ki9)))

    def v10(self, p, GAP, S7P, X5P, R5P):
        """ Reversible, fast EQ reaction: S7P + GAP <-> R5P

        (?) Sedoheptulose 7-phosphate + (?) Glyceraldehyde 3-phosphate
        -- Transketolase (TK)-->
        Ribulose 5-phosphate + Xylulose 5-phosphate
        """
        return p.kRE*((GAP*S7P) - (X5P*R5P)/p.q10)

    def v11(self, p, R5P, RU5P):
        """ Reversible, fast EQ reaction: R5P <-> Ru5P

        (?) Ribulose 5-phosphate
        -- Ribose 5-hosphate isomerase (RPI)-->
        (?) Ribulose 5-phosphate
        """
        return p.kRE*((R5P) - (RU5P)/p.q11)

    def v12(self, p, X5P, RU5P):
        """ Reversible, fast EQ reaction: X5P <-> Ru5P

        (?) Xylulose 5-phosphate
        -- Ribulose 5-phosphate 3 epimerase (RPE)-->
        (?) Ribulose 5-phosphate
        """
        return p.kRE*((X5P) - (RU5P)/p.q12)

    def v13(self, p, RU5P, ATP, PGA, RUBP, P, ADP):
        """ Irreversible reaction: Ru5P + ATP -> RuBP + ADP

        (3) Ribulose 5-phosphate + (3) ATP
        -- Phosphoribulokinase (PRK)-->
        (3) Ribulose 1, 5-bisphosphate + (3) ADP
        """
        return ((p.V13*RU5P*ATP)/((RU5P+p.Km131*(1+(PGA/p.Ki131)+(RUBP/p.Ki132)+(P/p.Ki133)))
                * (ATP*(1+(ADP/p.Ki134))+p.Km132*(1+(ADP/p.Ki135)))))

    def v14(self, p, F6P, G6P):
        """ Reversible, fast EQ reaction: F6P <-> G6P
        F6P
        -- Glucose 6-Phosphate isomerase (GPI)-->
        G6P
        """
        return p.kRE*((F6P) - (G6P)/p.q14)

    def v15(self, p, G6P, G1P):
        """ Reversible, fast EQ reaction: G6P <-> G1P
        G6P
        -- Phosphoglucomutase (PGM)
        --> G1P
        """
        return p.kRE*((G6P) - (G1P)/p.q15)

    def v16(self, p, ADP, Phosphate_i):
        """ATP regeneration  via ATP Synthase"""
        return (p.V16*ADP*Phosphate_i)/((ADP+p.Km161)*(Phosphate_i+p.Km162))

    def vPGA_out(self, p, PGA, N):
        """PGA export into medium"""
        return (p.Vx*PGA)/(N*p.Kpga)

    def vGAP_out(self, p, GAP, N):
        """GAP export into medium"""
        return (p.Vx*GAP)/(N*p.Kgap)

    def vDHAP_out(self, p, DHAP, N):
        """DHAP export into medium"""
        return (p.Vx*DHAP)/(N*p.Kdhap)

    def vStarchProduction(self, p, G1P, ATP, ADP, P, PGA, F6P, FBP):
        """G1P -> Gn-1 ; Starch production"""
        return ((p.Vst*G1P*ATP)
                / ((G1P+p.Kmst1)*((1+(ADP/p.Kist))*(ATP+p.Kmst2)
                   + ((p.Kmst2*P)/(p.Kast1*PGA+p.Kast2*F6P+p.Kast3*FBP)))))


class ReactionsNADPH(Reactions):
    ###########################################################################
    # Algebraic
    ###########################################################################
    def NADP(self, p, NADPH):
        return p.CN-NADPH

    ###########################################################################
    # Reaction rates
    ###########################################################################

    def v1(self, p, RUBP, PGA, FBP, SBP, P, NADPH):
        """
        3 Ribulose-1, 5-bisphosphate + 3 CO2
        -- RuBisCO -->
        6 3-Phosphoglycerate

        RuBp + CO2 -> PGA
        """
        return (p.V1*RUBP)/(RUBP+p.Km1*(1+(PGA/p.Ki11)+(FBP/p.Ki12)+(SBP/p.Ki13)+(P/p.Ki14)+(NADPH/p.Ki15)))

    def v3(self, p, NADPH, BPGA, GAP, NADP, P_i):
        """ Reversible, fast EQ reaction: BPGA + NADPH <-> GAP + NADP

        6 1, 3-Bisphosphoglycerate + 6 NADPH + 6 H+
        -- Glyceraldehyde 3-phosphate dehydrogenase (GADPH)-->
        1 G3P + 5 Glyceraldehyde 3-phosphate

        Stroma pH is assumed to be constant
        """
        return p.kRE*((NADPH*BPGA*p.protonsStroma) - (GAP*NADP*P_i)/p.q3)

    def vNADPH(self, p, NADP):
        return (p.Vnadph*NADP)/(p.Kmnadph+NADP)
