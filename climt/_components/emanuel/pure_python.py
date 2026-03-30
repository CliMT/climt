# -*- coding: utf-8 -*-
import numpy as np


class EmanuelConvectionPython(object):
    def __init__(self, **kwargs):
        # Default parameters from Fortran code
        self.IPBL = 0
        self.MINORIG = 1
        self.ELCRIT = 0.0011
        self.TLCRIT = -55.0
        self.ENTP = 1.5
        self.SIGD = 0.05
        self.SIGS = 0.12
        self.OMTRAIN = 50.0
        self.OMTSNOW = 5.5
        self.COEFFR = 1.0
        self.COEFFS = 0.8
        self.CU = 0.7
        self.BETA = 10.0
        self.DTMAX = 0.9
        self.ALPHA = 0.1
        self.DAMP = 0.1
        self.CPD = 1005.7
        self.CPV = 1870.0
        self.CL = 2500.0
        self.RV = 461.5
        self.RD = 287.04
        self.LV0 = 2.501e6
        self.G = 9.8
        self.ROWL = 1000.0
        self.DELT0 = 300.0

        # Update with kwargs
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def _tlift(self, P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW):
        # ND, NL, NK, KK, ICB are 0-based indices
        CPVMCL = self.CL - self.CPV
        EPS = self.RD / self.RV
        EPSI = 1.0 / EPS

        AH0 = (
            (self.CPD * (1.0 - Q[NK]) + self.CL * Q[NK]) * T[NK]
            + Q[NK] * (self.LV0 - CPVMCL * (T[NK] - 273.15))
            + GZ[NK]
        )
        CPP = self.CPD * (1.0 - Q[NK]) + Q[NK] * self.CPV
        CPINV = 1.0 / CPP

        if KK == 1:
            for i in range(ICB):
                CLW[i] = 0.0
            for i in range(NK, ICB):
                TPK[i] = T[NK] - (GZ[i] - GZ[NK]) * CPINV
                TVP[i] = TPK[i] * (1.0 + Q[NK] * EPSI)

        NST = ICB
        NSB = ICB
        if KK == 2:
            NST = NL
            NSB = ICB + 1

        for i in range(NSB, NST + 1):
            TG = T[i]
            QG = QS[i]
            ALV = self.LV0 - CPVMCL * (T[i] - 273.15)
            for j in range(2):
                S = self.CPD + ALV * ALV * QG / (self.RV * T[i] * T[i])
                S = 1.0 / S
                AHG = (
                    self.CPD * TG
                    + (self.CL - self.CPD) * Q[NK] * T[i]
                    + ALV * QG
                    + GZ[i]
                )
                TG = TG + S * (AH0 - AHG)
                TG = max(TG, 35.0)
                TC = TG - 273.15
                DENOM = 243.5 + TC
                if TC >= 0.0:
                    ES = 6.112 * np.exp(17.67 * TC / DENOM)
                else:
                    ES = np.exp(23.33086 - 6111.72784 / TG + 0.15215 * np.log(TG))
                QG = EPS * ES / (P[i] - ES * (1.0 - EPS))

            TPK[i] = (
                AH0 - (self.CL - self.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG
            ) / self.CPD
            CLW[i] = Q[NK] - QG
            CLW[i] = max(0.0, CLW[i])
            RG = QG / (1.0 - Q[NK])
            TVP[i] = TPK[i] * (1.0 + RG * EPSI)

        return TVP, TPK, CLW

    def _convect(
        self,
        T_in,
        Q_in,
        QS_in,
        U_in,
        V_in,
        P_in,
        PH_in,
        ND,
        NL,
        NTRA,
        DELT,
        CBMF_in,
        TRA_in,
    ):
        T = T_in.copy()
        Q = Q_in.copy()
        QS = QS_in.copy()
        U = U_in.copy()
        V = V_in.copy()
        P = P_in.copy()
        PH = PH_in.copy()
        TRA = TRA_in.copy()
        CBMF = CBMF_in

        FT = np.zeros(ND)
        FQ = np.zeros(ND)
        FU = np.zeros(ND)
        FV = np.zeros(ND)
        FTRA = np.zeros((ND, max(1, NTRA)))

        CPVMCL = self.CL - self.CPV
        EPS = self.RD / self.RV
        EPSI = 1.0 / EPS
        GINV = 1.0 / self.G
        DELTI = 1.0 / DELT

        TH = np.zeros(NL + 1)
        for i in range(NL + 1):
            RDCP = (self.RD * (1.0 - Q[i]) + Q[i] * self.RV) / (
                self.CPD * (1.0 - Q[i]) + Q[i] * self.CPV
            )
            TH[i] = T[i] * (1000.0 / P[i]) ** RDCP

        PRECIP = 0.0
        WD = 0.0
        TPRIME = 0.0
        QPRIME = 0.0
        IFLAG = 0
        OUTCAPE = 0.0

        if self.IPBL != 0:
            JC = 0
            for i in range(NL - 1, -1, -1):
                JN = -1
                SUM1 = TH[i] * (1.0 + Q[i] * EPSI - Q[i])
                for j in range(i + 1, NL):
                    SUM1 += TH[j] * (1.0 + Q[j] * EPSI - Q[j])
                    THBAR = SUM1 / float(j + 1 - i)
                    if (TH[j] * (1.0 + Q[j] * EPSI - Q[j])) < THBAR:
                        JN = j
                if i == 0:
                    JN = max(JN, 1)
                if JN == -1:
                    continue

                while True:
                    AHM = 0.0
                    RM = 0.0
                    UM = 0.0
                    VM = 0.0
                    TRATM = np.zeros(max(1, NTRA))
                    for j in range(i, JN + 1):
                        dp = PH[j] - PH[j + 1]
                        AHM += (self.CPD * (1.0 - Q[j]) + Q[j] * self.CPV) * T[j] * dp
                        RM += Q[j] * dp
                        UM += U[j] * dp
                        VM += V[j] * dp
                        for k in range(NTRA):
                            TRATM[k] += TRA[j, k] * dp

                    DPHINV = 1.0 / (PH[i] - PH[JN + 1])
                    RM *= DPHINV
                    UM *= DPHINV
                    VM *= DPHINV
                    for k in range(NTRA):
                        TRATM[k] *= DPHINV

                    A2 = 0.0
                    TOLD_slice = T[i : JN + 1].copy()
                    for j in range(i, JN + 1):
                        Q[j] = RM
                        U[j] = UM
                        V[j] = VM
                        for k in range(NTRA):
                            TRA[j, k] = TRATM[k]
                        RDCP = (self.RD * (1.0 - Q[j]) + Q[j] * self.RV) / (
                            self.CPD * (1.0 - Q[j]) + Q[j] * self.CPV
                        )
                        X = (0.001 * P[j]) ** RDCP
                        T[j] = X
                        A2 += (
                            (self.CPD * (1.0 - Q[j]) + Q[j] * self.CPV)
                            * X
                            * (PH[j] - PH[j + 1])
                        )

                    for idx, j in enumerate(range(i, JN + 1)):
                        TH_val = AHM / A2
                        T[j] = T[j] * TH_val
                        TH[j] = TH_val
                        TOLD_j = TOLD_slice[idx]
                        TC = TOLD_j - 273.15
                        ALV = self.LV0 - CPVMCL * TC
                        QS[j] = QS[j] + QS[j] * (1.0 + QS[j] * (EPSI - 1.0)) * ALV * (
                            T[j] - TOLD_j
                        ) / (self.RV * TOLD_j * TOLD_j)

                    if JN < NL - 1:
                        if (TH[JN + 1] * (1.0 + Q[JN + 1] * EPSI - Q[JN + 1])) < (
                            TH[JN] * (1.0 + Q[JN] * EPSI - Q[JN])
                        ):
                            JN += 1
                            continue
                    if i == 0:
                        JC = JN + 1
                    break

            if JC > 0:
                for j in range(JC):
                    if QS[j] < Q[j]:
                        ALV = self.LV0 - CPVMCL * (T[j] - 273.15)
                        denom_qs = (
                            self.CPD * (1.0 - Q[j])
                            + self.CL * Q[j]
                            + QS[j]
                            * (self.CPV - self.CL + ALV * ALV / (self.RV * T[j] * T[j]))
                        )
                        TNEW = T[j] + ALV * (Q[j] - QS[j]) / denom_qs
                        ALVNEW = self.LV0 - CPVMCL * (TNEW - 273.15)
                        QNEW = (
                            ALV * Q[j]
                            - (TNEW - T[j]) * (self.CPD * (1.0 - Q[j]) + self.CL * Q[j])
                        ) / ALVNEW
                        PRECIP += (
                            24.0
                            * 3600.0
                            * 1.0e5
                            * (PH[j] - PH[j + 1])
                            * (Q[j] - QNEW)
                            / (self.G * DELT * self.ROWL)
                        )
                        T[j] = TNEW
                        Q[j] = QNEW
                        QS[j] = QNEW

        GZ = np.zeros(ND + 1)
        CPN = np.zeros(ND + 1)
        H = np.zeros(ND + 1)
        LV = np.zeros(ND + 1)
        HM = np.zeros(ND + 1)
        TV = np.zeros(ND + 1)

        GZ[0] = 0.0
        CPN[0] = self.CPD * (1.0 - Q[0]) + Q[0] * self.CPV
        H[0] = T[0] * CPN[0]
        LV[0] = self.LV0 - CPVMCL * (T[0] - 273.15)
        HM[0] = LV[0] * Q[0]
        TV[0] = T[0] * (1.0 + Q[0] * EPSI - Q[0])

        AHMIN = 1.0e12
        IHMIN = NL

        for i in range(1, NL + 1):
            TVX = T[i] * (1.0 + Q[i] * EPSI - Q[i])
            TVY = T[i - 1] * (1.0 + Q[i - 1] * EPSI - Q[i - 1])
            GZ[i] = GZ[i - 1] + 0.5 * self.RD * (TVX + TVY) * (P[i - 1] - P[i]) / PH[i]
            CPN[i] = self.CPD * (1.0 - Q[i]) + self.CPV * Q[i]
            H[i] = T[i] * CPN[i] + GZ[i]
            LV[i] = self.LV0 - CPVMCL * (T[i] - 273.15)
            HM[i] = (
                (self.CPD * (1.0 - Q[i]) + self.CL * Q[i]) * (T[i] - T[0])
                + LV[i] * Q[i]
                + GZ[i]
            )
            TV[i] = T[i] * (1.0 + Q[i] * EPSI - Q[i])
            if (i + 1) >= self.MINORIG and HM[i] < AHMIN and HM[i] < HM[i - 1]:
                AHMIN = HM[i]
                IHMIN = i

        IHMIN = min(IHMIN, NL - 1)
        AHMAX = 0.0
        NK = 0
        for i in range(self.MINORIG - 1, IHMIN + 1):
            if HM[i] > AHMAX:
                NK = i
                AHMAX = HM[i]

        if T[NK] < 250.0 or Q[NK] <= 0.0 or IHMIN == (NL - 1):
            IFLAG = 0
            CBMF = 0.0
            return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

        RH = Q[NK] / QS[NK]
        CHI = T[NK] / (1669.0 - 122.0 * RH - T[NK])
        PLCL = P[NK] * (RH**CHI)

        if PLCL < 200.0 or PLCL >= 2000.0:
            IFLAG = 2
            CBMF = 0.0
            return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

        ICB = NL - 1
        for i in range(NK + 1, NL + 1):
            if P[i] < PLCL:
                ICB = min(ICB, i)

        if ICB >= (NL - 1):
            IFLAG = 3
            CBMF = 0.0
            return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

        TVP = np.zeros(ND)
        TP = np.zeros(ND)
        CLW = np.zeros(ND)
        TVP, TP, CLW = self._tlift(P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW)
        for i in range(NK, ICB + 1):
            TVP[i] -= TP[i] * Q[NK]

        if CBMF == 0.0 and TVP[ICB] <= (TV[ICB] - self.DTMAX):
            IFLAG = 0
            return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

        if IFLAG != 4:
            IFLAG = 1

        TVP, TP, CLW = self._tlift(P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW)

        EP = np.zeros(ND)
        SIGP = np.zeros(ND)
        for i in range(NK + 1):
            EP[i] = 0.0
            SIGP[i] = self.SIGS
        for i in range(NK + 1, NL + 1):
            TCA = TP[i] - 273.15
            ELACRIT = (
                self.ELCRIT if TCA >= 0.0 else self.ELCRIT * (1.0 - TCA / self.TLCRIT)
            )
            ELACRIT = max(ELACRIT, 0.0)
            EPMAX = 0.999
            EP[i] = EPMAX * (1.0 - ELACRIT / max(CLW[i], 1.0e-8))
            EP[i] = max(min(EP[i], EPMAX), 0.0)
            SIGP[i] = self.SIGS

        for i in range(ICB + 1, NL + 1):
            TVP[i] -= TP[i] * Q[NK]

        HP = H.copy()
        NENT = np.zeros(ND + 1, dtype=np.int32)
        WATER = np.zeros(ND + 1)
        EVAP = np.zeros(ND + 1)
        WT = np.full(ND + 1, self.OMTSNOW)
        MP = np.zeros(ND + 1)
        M = np.zeros(ND + 1)
        LVCP = np.zeros(ND + 1)
        LVCP[: NL + 1] = LV[: NL + 1] / CPN[: NL + 1]

        QENT = np.zeros((ND + 1, ND + 1))
        ELIJ = np.zeros((ND + 1, ND + 1))
        MENT = np.zeros((ND + 1, ND + 1))
        SIJ = np.zeros((ND + 1, ND + 1))
        UENT = np.zeros((ND + 1, ND + 1))
        VENT = np.zeros((ND + 1, ND + 1))
        TRAENT = np.zeros((ND + 1, ND + 1, max(1, NTRA)))

        for i in range(NL + 1):
            for j in range(NL + 1):
                QENT[i, j] = Q[j]
                UENT[i, j] = U[j]
                VENT[i, j] = V[j]
                for k in range(NTRA):
                    TRAENT[i, j, k] = TRA[j, k]

        QP = np.zeros(ND + 1)
        UP = np.zeros(ND + 1)
        VP = np.zeros(ND + 1)
        TRAP = np.zeros((ND + 1, max(1, NTRA)))

        QP[0] = Q[0]
        UP[0] = U[0]
        VP[0] = V[0]
        for i in range(NTRA):
            TRAP[0, i] = TRA[0, i]
        for i in range(1, NL + 1):
            QP[i] = Q[i - 1]
            UP[i] = U[i - 1]
            VP[i] = V[i - 1]
            for j in range(NTRA):
                TRAP[i, j] = TRA[i - 1, j]

        CAPE = 0.0
        CAPEM = 0.0
        INB = ICB + 1
        INB1 = INB
        BYP = 0.0
        for i in range(ICB + 1, NL):
            BY = (TVP[i] - TV[i]) * (PH[i] - PH[i + 1]) / P[i]
            CAPE += BY
            if BY >= 0.0:
                INB1 = i + 1
            if CAPE > 0.0:
                INB = i + 1
                BYP = (TVP[i + 1] - TV[i + 1]) * (PH[i + 1] - PH[i + 2]) / P[i + 1]
                CAPEM = CAPE
        INB = max(INB, INB1)
        CAPE = CAPEM + BYP
        DEFRAC = max(CAPEM - CAPE, 0.001)
        FRAC = min(max(-CAPE / DEFRAC, 0.0), 1.0)
        OUTCAPE = CAPE

        for i in range(ICB, INB + 1):
            HP[i] = H[NK] + (LV[i] + (self.CPD - self.CPV) * T[i]) * EP[i] * CLW[i]

        TVPPLCL = TVP[ICB - 1] - self.RD * TVP[ICB - 1] * (P[ICB - 1] - PLCL) / (
            CPN[ICB - 1] * P[ICB - 1]
        )
        TVAPLCL = TV[ICB] + (TVP[ICB] - TVP[ICB + 1]) * (PLCL - P[ICB]) / (
            P[ICB] - P[ICB + 1]
        )
        DTPBL = 0.0
        for i in range(NK, ICB):
            DTPBL += (TVP[i] - TV[i]) * (PH[i] - PH[i + 1])
        DTPBL /= PH[NK] - PH[ICB]
        DTMA = TVPPLCL - TVAPLCL + self.DTMAX + DTPBL

        CBMFOLD = CBMF
        DAMPS = self.DAMP * DELT / self.DELT0
        CBMF = (1.0 - DAMPS) * CBMF + 0.1 * self.ALPHA * DTMA
        CBMF = max(CBMF, 0.0)

        if CBMF == 0.0 and CBMFOLD == 0.0:
            return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

        M[ICB] = 0.0
        DBOSUM = 0.0
        for i in range(ICB + 1, INB + 1):
            k = min(i, INB1)
            DBO = abs(TV[k] - TVP[k]) + self.ENTP * 0.02 * (PH[k] - PH[k + 1])
            DBOSUM += DBO
            M[i] = CBMF * DBO
        if DBOSUM > 0:
            for i in range(ICB + 1, INB + 1):
                M[i] /= DBOSUM

        for i in range(ICB + 1, INB + 1):
            QTI = Q[NK] - EP[i] * CLW[i]
            for j in range(ICB, INB + 1):
                BF2 = 1.0 + LV[j] * LV[j] * QS[j] / (self.RV * T[j] * T[j] * self.CPD)
                ANUM = H[j] - HP[i] + (self.CPV - self.CPD) * T[j] * (QTI - Q[j])
                DENOM = H[i] - HP[i] + (self.CPD - self.CPV) * (Q[i] - QTI) * T[j]
                DEI = DENOM
                if abs(DEI) < 0.01:
                    DEI = 0.01
                SIJ[i, j] = ANUM / DEI
                SIJ[i, i] = 1.0
                ALTEM = (SIJ[i, j] * Q[i] + (1.0 - SIJ[i, j]) * QTI - QS[j]) / BF2
                CWAT = CLW[j] * (1.0 - EP[j])
                if (SIJ[i, j] < 0.0 or SIJ[i, j] > 1.0 or ALTEM > CWAT) and j > i:
                    ANUM -= LV[j] * (QTI - QS[j] - CWAT * BF2)
                    DENOM += LV[j] * (Q[i] - QTI)
                    if abs(DENOM) < 0.01:
                        DENOM = 0.01
                    SIJ[i, j] = ANUM / DENOM
                    ALTEM = (
                        SIJ[i, j] * Q[i]
                        + (1.0 - SIJ[i, j]) * QTI
                        - QS[j]
                        - (BF2 - 1.0) * CWAT
                    )

                if 0.0 < SIJ[i, j] < 0.9:
                    QENT[i, j] = SIJ[i, j] * Q[i] + (1.0 - SIJ[i, j]) * QTI
                    UENT[i, j] = SIJ[i, j] * U[i] + (1.0 - SIJ[i, j]) * U[NK]
                    VENT[i, j] = SIJ[i, j] * V[i] + (1.0 - SIJ[i, j]) * V[NK]
                    for k in range(NTRA):
                        TRAENT[i, j, k] = (
                            SIJ[i, j] * TRA[i, k] + (1.0 - SIJ[i, j]) * TRA[NK, k]
                        )
                    ELIJ[i, j] = max(0.0, ALTEM)
                    MENT[i, j] = M[i] / (1.0 - SIJ[i, j])
                    NENT[i] += 1
                SIJ[i, j] = min(max(SIJ[i, j], 0.0), 1.0)
            if NENT[i] == 0:
                MENT[i, i] = M[i]
                QENT[i, i] = Q[NK] - EP[i] * CLW[i]
                UENT[i, i] = U[NK]
                VENT[i, i] = V[NK]
                for k in range(NTRA):
                    TRAENT[i, i, k] = TRA[NK, k]
                ELIJ[i, i] = CLW[i]
                SIJ[i, i] = 1.0
        SIJ[INB, INB] = 1.0

        for i in range(ICB + 1, INB + 1):
            if NENT[i] != 0:
                QP1 = Q[NK] - EP[i] * CLW[i]
                ANUM = H[i] - HP[i] - LV[i] * (QP1 - QS[i])
                DENOM = H[i] - HP[i] + LV[i] * (Q[i] - QP1)
                if abs(DENOM) < 0.01:
                    DENOM = 0.01
                SCRIT = ANUM / DENOM
                ALT = QP1 - QS[i] + SCRIT * (Q[i] - QP1)
                if ALT < 0.0:
                    SCRIT = 1.0
                SCRIT = max(SCRIT, 0.0)
                ASIJ = 0.0
                SMIN = 1.0
                for j in range(ICB, INB + 1):
                    if 0.0 < SIJ[i, j] < 0.9:
                        if j > i:
                            SMID = min(SIJ[i, j], SCRIT)
                            SJMAX = SMID
                            SJMIN = SMID
                            if SMID < SMIN and SIJ[i, j + 1] < SMID:
                                SMIN = SMID
                                SJMAX = min(SIJ[i, j + 1], SIJ[i, j], SCRIT)
                                SJMIN = min(max(SIJ[i, j - 1], SIJ[i, j]), SCRIT)
                        else:
                            SJMAX = max(SIJ[i, j + 1], SCRIT)
                            SMID = max(SIJ[i, j], SCRIT)
                            SJMIN = max(SIJ[i, j - 1] if j > 0 else 0.0, SCRIT)
                        ASIJ += (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                            PH[j] - PH[j + 1]
                        )
                        MENT[i, j] *= (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                            PH[j] - PH[j + 1]
                        )
                ASIJ = max(1.0e-21, ASIJ)
                for j in range(ICB, INB + 1):
                    MENT[i, j] /= ASIJ
                if sum(MENT[i, ICB : INB + 1]) < 1.0e-18:
                    NENT[i] = 0
                    MENT[i, i] = M[i]
                    QENT[i, i] = Q[NK] - EP[i] * CLW[i]
                    UENT[i, i] = U[NK]
                    VENT[i, i] = V[NK]
                    for k in range(NTRA):
                        TRAENT[i, i, k] = TRA[NK, k]
                    ELIJ[i, i] = CLW[i]
                    SIJ[i, i] = 1.0

        if EP[INB] >= 0.0001:
            JTT = 0
            for i in range(INB, -1, -1):
                WDTRAIN = self.G * EP[i] * M[i] * CLW[i]
                if i > 0:
                    for j in range(i):
                        WDTRAIN += (
                            self.G
                            * max(0.0, ELIJ[j, i] - (1.0 - EP[i]) * CLW[i])
                            * MENT[j, i]
                        )
                COEFF = self.COEFFR if T[i] > 273.0 else self.COEFFS
                WT[i] = self.OMTRAIN if T[i] > 273.0 else self.OMTSNOW
                AFAC = max(
                    COEFF
                    * PH[i]
                    * (QS[i] - 0.5 * (Q[i] + QP[i + 1]))
                    / (1.0e4 + 2.0e3 * PH[i] * QS[i]),
                    0.0,
                )
                SIGT = min(max(SIGP[i], 0.0), 1.0)
                B6 = 100.0 * (PH[i] - PH[i + 1]) * SIGT * AFAC / WT[i]
                C6 = (WATER[i + 1] * WT[i + 1] + WDTRAIN / self.SIGD) / WT[i]
                REVAP = 0.5 * (-B6 + np.sqrt(B6 * B6 + 4.0 * C6))
                EVAP[i] = SIGT * AFAC * REVAP
                WATER[i] = REVAP * REVAP
                if i > 0:
                    DHDP = max((H[i] - H[i - 1]) / (P[i - 1] - P[i]), 10.0)
                    MP[i] = 100.0 * GINV * LV[i] * self.SIGD * EVAP[i] / DHDP
                    FAC = 20.0 / (PH[i - 1] - PH[i])
                    MP[i] = (FAC * MP[i + 1] + MP[i]) / (1.0 + FAC)
                    if P[i] > (0.949 * P[0]):
                        JTT = max(JTT, i)
                        MP[i] = MP[JTT] * (P[0] - P[i]) / (P[0] - P[JTT])
                if i != INB:
                    QSTM = QS[0] if i == 0 else QS[i - 1]
                    if MP[i] > MP[i + 1]:
                        RAT = MP[i + 1] / MP[i]
                        QP[i] = (
                            QP[i + 1] * RAT
                            + Q[i] * (1.0 - RAT)
                            + 100.0
                            * GINV
                            * self.SIGD
                            * (PH[i] - PH[i + 1])
                            * (EVAP[i] / MP[i])
                        )
                        UP[i] = UP[i + 1] * RAT + U[i] * (1.0 - RAT)
                        VP[i] = VP[i + 1] * RAT + V[i] * (1.0 - RAT)
                        for k in range(NTRA):
                            TRAP[i, k] = TRAP[i + 1, k] * RAT + TRAP[i, k] * (1.0 - RAT)
                    elif MP[i + 1] > 0.0:
                        QP[i] = (
                            GZ[i + 1]
                            - GZ[i]
                            + QP[i + 1] * (LV[i + 1] + T[i + 1] * (self.CL - self.CPD))
                            + self.CPD * (T[i + 1] - T[i])
                        ) / (LV[i] + T[i] * (self.CL - self.CPD))
                        UP[i] = UP[i + 1]
                        VP[i] = VP[i + 1]
                        for k in range(NTRA):
                            TRAP[i, k] = TRAP[i + 1, k]
                    QP[i] = max(min(QP[i], QSTM), 0.0)
            PRECIP += (
                WT[0] * self.SIGD * WATER[0] * 3600.0 * 24000.0 / (self.ROWL * self.G)
            )

        WD = self.BETA * abs(MP[ICB]) * 0.01 * self.RD * T[ICB] / (self.SIGD * P[ICB])
        QPRIME = 0.5 * (QP[0] - Q[0])
        TPRIME = self.LV0 * QPRIME / self.CPD
        DPINV = 0.01 / (PH[0] - PH[1])
        AM = 0.0
        if NK == 0:
            AM = sum(M[1 : INB + 1])
        if (2.0 * self.G * DPINV * AM) >= DELTI:
            IFLAG = 4
        FT[0] += (
            self.G * DPINV * AM * (T[1] - T[0] + (GZ[1] - GZ[0]) / CPN[0])
            - LVCP[0] * self.SIGD * EVAP[0]
            + self.SIGD
            * WT[1]
            * (self.CL - self.CPD)
            * WATER[1]
            * (T[1] - T[0])
            * DPINV
            / CPN[0]
        )
        FQ[0] += (
            self.G * MP[1] * (QP[1] - Q[0]) * DPINV
            + self.SIGD * EVAP[0]
            + self.G * AM * (Q[1] - Q[0]) * DPINV
        )
        FU[0] += self.G * DPINV * (MP[1] * (UP[1] - U[0]) + AM * (U[1] - U[0]))
        FV[0] += self.G * DPINV * (MP[1] * (VP[1] - V[0]) + AM * (V[1] - V[0]))
        for j in range(NTRA):
            FTRA[0, j] += (
                self.G
                * DPINV
                * (MP[1] * (TRAP[1, j] - TRA[0, j]) + AM * (TRA[1, j] - TRA[0, j]))
            )
        for j in range(1, INB + 1):
            FQ[0] += self.G * DPINV * MENT[j, 0] * (QENT[j, 0] - Q[0])
            FU[0] += self.G * DPINV * MENT[j, 0] * (UENT[j, 0] - U[0])
            FV[0] += self.G * DPINV * MENT[j, 0] * (VENT[j, 0] - V[0])
            for k in range(NTRA):
                FTRA[0, k] += (
                    self.G * DPINV * MENT[j, 0] * (TRAENT[j, 0, k] - TRA[0, k])
                )

        for i in range(1, INB + 1):
            DPINV = 0.01 / (PH[i] - PH[i + 1])
            CPINV = 1.0 / CPN[i]
            AMP1 = sum(M[i + 1 : INB + 2]) if i >= NK else 0.0
            for k in range(i + 1):
                AMP1 += sum(MENT[k, i + 1 : INB + 2])
            if (2.0 * self.G * DPINV * AMP1) >= DELTI:
                IFLAG = 4
            AD = 0.0
            for k in range(i):
                AD += sum(MENT[i : INB + 1, k])
            FT[i] += (
                self.G
                * DPINV
                * (
                    AMP1 * (T[i + 1] - T[i] + (GZ[i + 1] - GZ[i]) * CPINV)
                    - AD * (T[i] - T[i - 1] + (GZ[i] - GZ[i - 1]) * CPINV)
                )
                - self.SIGD * LVCP[i] * EVAP[i]
            )
            FT[i] += (
                self.G
                * DPINV
                * MENT[i, i]
                * (HP[i] - H[i] + T[i] * (self.CPV - self.CPD) * (Q[i] - QENT[i, i]))
                * CPINV
            )
            FT[i] += (
                self.SIGD
                * WT[i + 1]
                * (self.CL - self.CPD)
                * WATER[i + 1]
                * (T[i + 1] - T[i])
                * DPINV
                * CPINV
            )
            FQ[i] += (
                self.G * DPINV * (AMP1 * (Q[i + 1] - Q[i]) - AD * (Q[i] - Q[i - 1]))
            )
            FU[i] += (
                self.G * DPINV * (AMP1 * (U[i + 1] - U[i]) - AD * (U[i] - U[i - 1]))
            )
            FV[i] += (
                self.G * DPINV * (AMP1 * (V[i + 1] - V[i]) - AD * (V[i] - V[i - 1]))
            )
            for k in range(NTRA):
                FTRA[i, k] += (
                    self.G
                    * DPINV
                    * (
                        AMP1 * (TRA[i + 1, k] - TRA[i, k])
                        - AD * (TRA[i, k] - TRA[i - 1, k])
                    )
                )
            for k in range(i):
                AWAT = max(ELIJ[k, i] - (1.0 - EP[i]) * CLW[i], 0.0)
                FQ[i] += self.G * DPINV * MENT[k, i] * (QENT[k, i] - AWAT - Q[i])
                FU[i] += self.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])
                FV[i] += self.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i])
                for j in range(NTRA):
                    FTRA[i, j] += (
                        self.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j])
                    )
            for k in range(i, INB + 1):
                FQ[i] += self.G * DPINV * MENT[k, i] * (QENT[k, i] - Q[i])
                FU[i] += self.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])
                FV[i] += self.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i])
                for j in range(NTRA):
                    FTRA[i, j] += (
                        self.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j])
                    )
            FQ[i] += (
                self.SIGD * EVAP[i]
                + self.G
                * (MP[i + 1] * (QP[i + 1] - Q[i]) - MP[i] * (QP[i] - Q[i - 1]))
                * DPINV
            )
            FU[i] += (
                self.G
                * (MP[i + 1] * (UP[i + 1] - U[i]) - MP[i] * (UP[i] - U[i - 1]))
                * DPINV
            )
            FV[i] += (
                self.G
                * (MP[i + 1] * (VP[i + 1] - V[i]) - MP[i] * (VP[i] - V[i - 1]))
                * DPINV
            )
            for j in range(NTRA):
                FTRA[i, j] += (
                    self.G
                    * DPINV
                    * (
                        MP[i + 1] * (TRAP[i + 1, j] - TRA[i, j])
                        - MP[i] * (TRAP[i, j] - TRA[i - 1, j])
                    )
                )

        FQOLD = FQ[INB]
        FQ[INB] *= 1.0 - FRAC
        FQ[INB - 1] += (
            FRAC
            * FQOLD
            * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
            * LV[INB]
            / LV[INB - 1]
        )
        FTOLD = FT[INB]
        FT[INB] *= 1.0 - FRAC
        FT[INB - 1] += (
            FRAC
            * FTOLD
            * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
            * CPN[INB]
            / CPN[INB - 1]
        )
        FUOLD = FU[INB]
        FU[INB] *= 1.0 - FRAC
        FU[INB - 1] += (
            FRAC * FUOLD * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
        )
        FVOLD = FV[INB]
        FV[INB] *= 1.0 - FRAC
        FV[INB - 1] += (
            FRAC * FVOLD * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
        )
        for k in range(NTRA):
            FTRAOLD = FTRA[INB, k]
            FTRA[INB, k] *= 1.0 - FRAC
            FTRA[INB - 1, k] += (
                FRAC * FTRAOLD * (PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB])
            )

        ENTS = sum(
            (CPN[: INB + 1] * FT[: INB + 1] + LV[: INB + 1] * FQ[: INB + 1])
            * (PH[: INB + 1] - PH[1 : INB + 2])
        ) / (PH[0] - PH[INB + 1])
        UAV = sum(FU[: INB + 1] * (PH[: INB + 1] - PH[1 : INB + 2])) / (
            PH[0] - PH[INB + 1]
        )
        VAV = sum(FV[: INB + 1] * (PH[: INB + 1] - PH[1 : INB + 2])) / (
            PH[0] - PH[INB + 1]
        )
        for i in range(INB + 1):
            FT[i] -= ENTS / CPN[i]
            FU[i] = (1.0 - self.CU) * (FU[i] - UAV)
            FV[i] = (1.0 - self.CU) * (FV[i] - VAV)
        for k in range(NTRA):
            TRAAV = sum(FTRA[: INB + 1, k] * (PH[: INB + 1] - PH[1 : INB + 2])) / (
                PH[0] - PH[INB + 1]
            )
            FTRA[: INB + 1, k] -= TRAAV

        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    def array_call(self, state, timestep):
        t = state["air_temperature"]
        q = state["specific_humidity"]
        u = state["eastward_wind"]
        v = state["northward_wind"]
        p = state["air_pressure"]
        ph = state["air_pressure_on_interface_levels"]

        nlev, ncol = t.shape

        try:
            from climt._core import bolton_q_sat
        except ImportError:

            def bolton_q_sat(T, p, Rd, Rv):
                tc = T - 273.15
                es = 611.2 * np.exp(17.67 * tc / (tc + 243.5))
                return (Rd / Rv) * es / (p - es * (1 - Rd / Rv))

        qs = bolton_q_sat(t, p * 100, self.RD, self.RV)

        ft = np.zeros_like(t)
        fq = np.zeros_like(q)
        fu = np.zeros_like(u)
        fv = np.zeros_like(v)
        precip = np.zeros(ncol)
        wd = np.zeros(ncol)
        tprime = np.zeros(ncol)
        qprime = np.zeros(ncol)
        outcape = np.zeros(ncol)
        iflag = np.zeros(ncol, dtype=np.int32)

        cbmf = state.get("cloud_base_mass_flux", np.zeros(ncol)).copy()
        ntra = 0
        tra = np.zeros((nlev, 1))
        delt = timestep.total_seconds()

        for col in range(ncol):
            (
                ft[:, col],
                fq[:, col],
                fu[:, col],
                fv[:, col],
                precip[col],
                wd[col],
                tprime[col],
                qprime[col],
                cbmf[col],
                outcape[col],
                iflag[col],
            ) = self._convect(
                t[:, col],
                q[:, col],
                qs[:, col],
                u[:, col],
                v[:, col],
                p[:, col],
                ph[:, col],
                nlev,
                nlev - 3,
                ntra,
                delt,
                cbmf[col],
                tra,
            )

        tendencies = {
            "air_temperature": ft,
            "specific_humidity": fq,
            "eastward_wind": fu,
            "northward_wind": fv,
        }
        diagnostics = {
            "convective_state": iflag,
            "convective_precipitation_rate": precip,
            "convective_downdraft_velocity_scale": wd,
            "convective_downdraft_temperature_scale": tprime,
            "convective_downdraft_specific_humidity_scale": qprime,
            "cloud_base_mass_flux": cbmf,
            "atmosphere_convective_available_potential_energy": outcape,
        }
        return tendencies, diagnostics
