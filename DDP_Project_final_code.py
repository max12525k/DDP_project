# Code of the program developed for this thesis

# Modules Used in this program

import math

import numpy as np
import pandas as pd
import rdflib
import scipy
import scipy.stats
from owlready2 import *
from sklearn import linear_model

# In[3]: Reading owl file


onto = get_ontology("file:///home/max/DDP_project.owl")
onto.load()


# In[4]: reading loading files


data_Shaft = pd.read_excel("Loading.xlsx")
data_Gear = pd.read_excel("Loading_Gear.xlsx")

# In[6]: compute reliablity


def reliability(TTFA, Beta, TTFT):
    rel = math.exp(-((TTFA / TTFT) ** Beta))
    return rel


# calculate moment
def Shaft_Moment(Force, Diameter, Length):  # calculate moment
    mom = Force * Length / 2
    return mom


# In[7]: does the linear regression
def linearreg(df, YA, Tn):
    X = np.array([])
    Y = np.array([])
    for i in range(len(df)):
        X = np.append(X, math.log(df[i]))
        Y = np.append(Y, math.log(-math.log(YA[i])))
    X = X.reshape(-1, 1)
    lm = linear_model.LinearRegression()
    model = lm.fit(X, Y)
    a = float(lm.coef_)
    b = float(lm.intercept_)
    neta = math.exp(-b / a)
    print(
        "Linear regression model for reliability and time to failure has yielded with Shape Parameter Neta = ",
        neta,
        "and Scale Parameter Beta = ",
        a,
    )
    t = math.exp(((math.log(-math.log(Tn))) + a * math.log(neta)) / a)
    t = int(t / 60)
    return t


# computes von mises stress
def Shaft_SigmaRev(Force, Torque, Diameter, Length):
    Dia = 3.14 * Diameter * Diameter * Diameter
    mom = Shaft_Moment(Force, Diameter, Length)
    Ma = mom * 32 / Dia
    ta = Torque * 16 / Dia
    sigmarev = math.sqrt(Ma * Ma + 3 * ta * ta)
    return sigmarev


# In[8]:compute f
def Shaft_f(UltimateTensileStrength):
    f = 0.0
    const = (
        80 * 0.876
        + 100 * 0.842
        + 120 * 0.82
        + 0.805 * 140
        + 0.792 * 160
        + 0.781 * 180
        + 0.775 * 200
    ) / 7
    if UltimateTensileStrength <= 70:
        f = 0.9
    else:
        f = const / UltimateTensileStrength
    return f


# In[9]:
# compute moment reliability


def Shaft_Relmoment(Force, YeildStrength, Diameter, Length, StdYsratio, StdFratio):
    mom = Shaft_Moment(Force, Diameter, Length)
    Dia = 3.14 * Diameter * Diameter * Diameter
    mom = mom * 32 / Dia
    Ys = YeildStrength / 2
    Mean = Ys - mom
    stdYs = StdYsratio * Ys
    stdF = StdFratio * mom
    Std = math.sqrt(stdYs * stdYs + stdF * stdF)
    norm = scipy.stats.norm(Mean, Std)
    rel = 1 - norm.cdf(-Mean / Std)
    return rel


# In[10]:compute torsion reliability


def Shaft_Reltorsion(Torque, YeildStrength, Diameter, Length, StdYsratio, StdTratio):
    Dia = 3.14 * Diameter * Diameter * Diameter
    Tor = Torque * 16 / Dia
    Ys = YeildStrength / 2
    Mean = Ys - Tor
    stdYs = StdYsratio * Ys
    stdT = StdTratio * Tor
    Std = math.sqrt(stdYs * stdYs + stdT * stdT)
    norm = scipy.stats.norm(Mean, Std)
    rel = 1 - norm.cdf(-Mean / Std)
    return rel


# In[11]:


def Shaft_Ka(UltimateTensileStrengthMPA, ManufactoringProcess):
    a = 0
    b = 0
    if ManufactoringProcess == "Ground":
        a = 1.58
        b = -0.085

    elif ManufactoringProcess == "Colddrawn" or ManufactoringProcess == "Machined":
        a = 4.51
        b = -0.265
    elif ManufactoringProcess == "Hotrolled":
        a = 57.7
        b = -0.718
    else:
        a = 272
        b = -0.995
    ka = a * UltimateTensileStrengthMPA ** b

    return ka


# In[12]:


def Shaft_Kb(Diameter):
    kb = 0
    if Diameter <= 51:
        kb = 1.24 * (Diameter) ** (-0.107)
    else:
        kb = 1.51 * (Diameter) ** (-0.157)
    return kb


# In[13]:


def Shaft_Kd(Temp):
    kd = (
        0.975
        + 0.432 * Temp * (10 ** (-3))
        - 0.115 * (10 ** (-5)) * Temp ** 2
        + 0.104 * (10 ** (-8)) * Temp ** 3
        - 0.595 * (10 ** (-12)) * Temp ** 4
    )
    return kd


# In[14]:


def Shaft_Ke(Za):
    if Za == 99:
        a = 2.326
    elif Za == 99.9:
        a = 3.091
    elif Za == 99.99:
        a = 3.719
    elif Za == 99.999:
        a = 4.265
    else:
        a = 0

    Ke = 1 - 0.08 * a
    return Ke


# In[15]: compute modified endurance strength


def Shaft_ModEnd(
    EnduranceLimit, UltimateTensileStrength, ManufactoringProcess, Diameter, Temp, Za
):
    ka = Shaft_Ka(UltimateTensileStrength, ManufactoringProcess)
    kb = Shaft_Kb(Diameter)
    kc = 0.85
    kd = Shaft_Kd(Temp)
    ke = Shaft_Ke(Za)
    ModifiedEnduranceLimit = ka * kb * kc * kd * ke * EnduranceLimit
    return ModifiedEnduranceLimit


# In[16]:life of the shaft


def Shaft_RelLife(
    RPM,
    Force,
    Torque,
    Length,
    EnduranceLimit,
    UltimateTensileStrength,
    ManufactoringProcess,
    Diameter,
    Temp,
    Za,
):
    modend = Shaft_ModEnd(
        EnduranceLimit,
        UltimateTensileStrength,
        ManufactoringProcess,
        Diameter,
        Temp,
        Za,
    )
    fc = Shaft_f(UltimateTensileStrength)
    a = fc * UltimateTensileStrength * fc * UltimateTensileStrength / modend
    b = math.log10(modend / fc / UltimateTensileStrength) / 3
    sigrev = Shaft_SigmaRev(Force, Torque, Diameter, Length)
    cycle = (sigrev / a) ** (1 / b)
    durn = cycle / RPM
    return durn


# in[16]
# In[6]:pitch diameter


def Gear_PitchDia(DiametralPitch, NoofTeeth):
    pitchdia = float(NoofTeeth) / DiametralPitch
    return pitchdia


# In[7]:allowable bending stress
def Gear_ST(Material, Heattreatment, Bhardness):
    St = 0.0
    if Material == "AISI4140" and Heattreatment == "Nitrided":
        St = 82.3 * Bhardness + 12150
    elif Material == "AISI4340" and Heattreatment == "Nitrided":
        St = 108.6 * Bhardness + 15890
    elif Material == "2.5%chrome" and Heattreatment == "Nitrided":
        St = 105.2 * Bhardness + 9280
    else:
        St = 102 * Bhardness + 16400
    return St


# In[8]:allowable contact stress


def Gear_Sc(Material, Heattreatment, Bhardness):
    Sc = 0.0
    if Material == "AISI4140" and Heattreatment == "Nitrided":
        Sc = 150000
    elif Material == "AISI4340" and Heattreatment == "Nitrided":
        St = 163000
    elif Material == "2.5%chrome" and Heattreatment == "Nitrided":
        Sc = 155000
    else:
        Sc = 349 * Bhardness + 34300
    return Sc


# In[9]: compute gear velocity
def Gear_velocity(Revolution, DiametralPitchp, NoofTeethp):
    Dp = Gear_PitchDia(DiametralPitchp, NoofTeethp)
    Vel = 3.14 * Dp * Revolution / 12
    return Vel


# In[10]:dynamic factor
def Gear_Kv(Qualitynumber, Revolution, DiametralPitchp, NoofTeethp):
    Vel = Gear_velocity(Revolution, DiametralPitchp, NoofTeethp)
    B = 0.25 * (12 - Qualitynumber) ** (2 / 3)
    A = 5.0 + 56 * (1 - B)
    Kv = ((A + (Vel) ** (0.5)) / A) ** B
    return Kv


# In[11]:size factor
def Gear_Ks(FaceWidth, LewisFormFactor, DiametralPitchg):
    F = FaceWidth
    Y = LewisFormFactor
    Ks = 1.192 * (F * (Y) ** (0.5) / DiametralPitchg)
    if Ks > 1:
        return Ks
    else:
        return 1


# In[12]: load distribution factor


def Gear_Km(TeethType, FaceWidth, DiametralPitchg, NoofTeethg):
    Dg = Gear_PitchDia(DiametralPitchg, NoofTeethg)
    Cmc = 0.0
    Cpf = 0.0
    Cpm = 1.0
    Ce = 1.0
    Cma = 0.0
    # for commercial enclosed unit
    A = 0.127
    B = 0.0158
    C = -0.930 / 10000
    Cma = A + B * FaceWidth + C * FaceWidth ** 2
    if TeethType == "Uncrowned":
        Cmc = 1
    else:
        Cmc = 0.8
    if FaceWidth <= 1:
        Cpf = FaceWidth / (10 * Dg) - 0.025
    elif FaceWidth <= 17 and FaceWidth >= 1:
        Cpf = FaceWidth / (10 * Dg) - 0.0375 + 0.0125 * FaceWidth
    else:
        Cpf = (
            FaceWidth / (10 * Dg)
            + 0.1109
            + 0.0207 * FaceWidth
            - 0.000228 * FaceWidth ** 2
        )
    Km = 1 + Cmc * (Cpf * Cpm + Cma * Ce)
    return Km


# In[13]:rim thickness factor


def Gear_Kb(ToothHight, RimThickness):
    Mb = RimThickness / ToothHight
    if Mb < 1.2:
        Kb = 1.6 * math.log(2.242 / Mb)
    else:
        Kb = 1
    return Kb


# In[14]:reliability factor
def Gear_Kr(Reliability):
    if Reliability >= 0.5 and Reliability < 0.99:
        Kr = 0.658 - 0.0759 * math.log(1 - Reliability)
    else:
        Kr = 0.5 - 0.109 * math.log(1 - Reliability)
    return Kr


# In[15]:compute curvature sum
def Gear_Rhosum(
    PressureAngle,
    NoofTeethg,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
):
    Dp = Gear_PitchDia(DiametralPitchp, NoofTeethp)
    Dg = Gear_PitchDia(DiametralPitchg, NoofTeethg)
    Rp = Dp / 2
    Rg = Dg / 2
    Bcrg = Rg * math.cos(PressureAngle)
    Bcrp = Rp * math.cos(PressureAngle)
    Z = (
        math.sqrt(AddendumRadiusg ** 2 - Bcrg ** 2)
        + math.sqrt(AddendumRadiusp ** 2 - Bcrp ** 2)
        - (Rp + Rg) * math.sin(PressureAngle)
    )
    Pb = 2 * 3.14 * Bcrg / NoofTeethg
    BetaL = (Z - Pb) / Bcrg
    BetaH = (2 * Pb - Z) / Bcrg
    Delta = (
        (Rp + Rg) * math.sin(PressureAngle)
        - math.sqrt(AddendumRadiusp ** 2 - Bcrp ** 2)
    ) / Bcrg
    Rocg = Bcrg * (Delta + BetaL)
    Rocp = (Rg + Rp) * math.sin(PressureAngle) - Rocg
    RocSum = (1 / Rocg) + (1 / Rocp)
    return RocSum


# In[16]:hardness ratio factor
def Gear_Ch(BHardnessg, BHardnessp, DiametralPitchg, DiametralPitchp):
    Mg = DiametralPitchg / DiametralPitchp
    Hbpg = BHardnessg / BHardnessp
    if Hbpg < 1.2:
        A = 0.0
    elif Hbpg > 1.7:
        A = 0.00698
    else:
        A = Hbpg * 8.98 * 10 ** (-3) - 8.29 * 10 ** (-3)
    Ch = 1.0 + A * (Mg - 1.0)
    return Ch


# In[17]:compute life of gear
def Gear_life(
    Load,
    PressureAngle,
    NoofTeethg,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
):
    Dp = Gear_PitchDia(DiametralPitchp, NoofTeethp)
    Dg = Gear_PitchDia(DiametralPitchg, NoofTeethg)
    Rp = Dp / 2
    Rg = Dg / 2
    Bcrg = Rg * math.cos(PressureAngle)
    Bcrp = Rp * math.cos(PressureAngle)
    Z = (
        math.sqrt(AddendumRadiusg ** 2 - Bcrg ** 2)
        + math.sqrt(AddendumRadiusp ** 2 - Bcrp ** 2)
        - (Rp + Rg) * math.sin(PressureAngle)
    )
    Pb = 2 * 3.14 * Bcrg / NoofTeethg
    BetaL = (Z - Pb) / Bcrg
    BetaH = (2 * Pb - Z) / Bcrg
    Delta = (
        (Rp + Rg) * math.sin(PressureAngle)
        - math.sqrt(AddendumRadiusp ** 2 - Bcrp ** 2)
    ) / Bcrg
    InvL = Bcrg * BetaH * (Delta + BetaL + BetaH / 2)
    Rocg = Bcrg * (Delta + BetaL)
    Rocp = (Rg + Rp) * math.sin(PressureAngle) - Rocg
    RocSum = (1 / Rocg) + (1 / Rocp)
    T10 = (
        3.72
        * 10 ** 18
        * Load ** (-4.3)
        * FaceWidth ** 3.9
        * RocSum ** (-5)
        * InvL ** (-0.4)
    )
    G10 = (NoofTeethg * (1 / T10) ** 2.5) ** (-1 / 2.5) * 10 ** 6
    return G10


# In[18]:compute stress cycle factor
def Gear_Yn(
    Load,
    PressureAngle,
    NoofTeethg,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
    Heattreatment,
    Bhardness,
):
    N = Gear_life(
        Load,
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
    )
    if Bhardness == 160 and N < 10 ** 6:
        Yn = 2.3194 * N ** (-0.0583)
    elif Bhardness == 250 and N < 10 ** 6:
        Yn = 4.9404 * N * (-0.1045)
    elif Bhardness == 400 and N < 10 ** 6:
        Yn = 9.4518 * N ** (-0.148)
    elif Heattreatment == "Nitrided" and N < 10 ** 6:
        Yn = 3.517 * N ** (-0.0817)
    else:
        Yn = 1.3558 * N ** (-0.0178)
    return Yn


# In[19]:stress cycle factor
def Gear_Zn(
    Load,
    PressureAngle,
    NoofTeethg,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
    Heattreatment,
):
    N = Gear_life(
        Load,
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
    )
    if Heattreatment == "Nitrided" and N < 10 ** 6:
        Zn = 1.249 * N ** (-0.0138)
    elif N < 10 ** 6:
        Zn = 2.466 * N ** (-0.056)
    else:
        Zn = 1.4488 * N ** (-0.023)
    return Zn


# In[20]:overload factor
def Gear_Ko(imptdriver, imptmachine):
    Ko = 0.0
    if imptdriver == "Uniform" and imptmachine == "Uniform":
        Ko = 1
    elif imptdriver == "Uniform" and imptmachine == "Moderate":
        Ko = 1.25
    elif imptdriver == "Uniform" and imptmachine == "Heavy":
        Ko = 1.75
    elif imptdriver == "Light" and imptmachine == "Uniform":
        Ko = 1.25
    elif imptdriver == "Light" and imptmachine == "Moderate":
        Ko = 1.50
    elif imptdriver == "Light" and imptmachine == "Heavy":
        Ko = 2.0
    elif imptdriver == "Medium" and imptmachine == "Uniform":
        Ko = 1.50
    elif imptdriver == "Medium" and imptmachine == "Moderate":
        Ko = 1.75
    elif imptdriver == "Medium" and imptmachine == "Heavy":
        Ko = 2.25
    else:
        print("Overload factor Not available")
    return Ko


# In[21]:Reliability bending
def Gear_RelB(
    TeethType,
    NoofTeethg,
    Load,
    StdStratio,
    StdLratio,
    BSGFactor,
    Qualitynumber,
    Bhardness,
    Reliability,
    Revolution,
    LewisFormFactor,
    RimThickness,
    ToothHight,
    imptdriver,
    imptmachine,
    PressureAngle,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
    Heattreatment,
    Material,
):
    ko = Gear_Ko(imptdriver, imptmachine)
    kv = Gear_Kv(Qualitynumber, Revolution, DiametralPitchp, NoofTeethp)
    ks = Gear_Ks(FaceWidth, LewisFormFactor, DiametralPitchg)
    km = Gear_Km(TeethType, FaceWidth, DiametralPitchg, NoofTeethg)
    kb = Gear_Kb(ToothHight, RimThickness)
    yn = Gear_Yn(
        Load,
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
        Heattreatment,
        Bhardness,
    )
    kt = 1
    kr = Gear_Kr(Reliability)
    St = Gear_ST(Material, Heattreatment, Bhardness)
    Sigmab = Load * ko * kv * ks * km * kb * DiametralPitchp / (FaceWidth * BSGFactor)
    Sta = St * yn / (kt * kr)
    Mean = Sta - Sigmab
    stdSta = StdStratio * Sta
    stdL = StdLratio * Sigmab
    Std = math.sqrt(stdSta * stdSta + stdL * stdL)
    norm = scipy.stats.norm(Mean, Std)
    rel = 1 - norm.cdf(-Mean / Std)
    return rel


# In[22]:reliability contact
def Gear_RelC(
    NoofTeethg,
    TeethType,
    Material,
    Load,
    ElasticCoefficient,
    StdScratio,
    StdLratio,
    BSGFactor,
    Qualitynumber,
    Bhardness,
    Reliability,
    Revolution,
    LewisFormFactor,
    RimThickness,
    ToothHight,
    imptdriver,
    imptmachine,
    PressureAngle,
    NoofTeethp,
    FaceWidth,
    AddendumRadiusg,
    AddendumRadiusp,
    DiametralPitchg,
    DiametralPitchp,
    Heattreatment,
):
    Bhardnessg = Bhardness
    Bhardnessp = Bhardness
    ko = Gear_Ko(imptdriver, imptmachine)
    kv = Gear_Kv(Qualitynumber, Revolution, DiametralPitchp, NoofTeethp)
    ks = Gear_Ks(FaceWidth, LewisFormFactor, DiametralPitchg)
    km = Gear_Km(TeethType, FaceWidth, DiametralPitchg, NoofTeethg)
    zn = Gear_Zn(
        Load,
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
        Heattreatment,
    )
    cf = 1
    kt = 1
    kr = Gear_Kr(Reliability)
    ch = Gear_Ch(Bhardnessg, Bhardnessp, DiametralPitchg, DiametralPitchp)
    rhosum = Gear_Rhosum(
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
    )
    Sc = Gear_Sc(Material, Heattreatment, Bhardness)
    Sigmac = ElasticCoefficient * (
        math.sqrt(abs(Load * ko * kv * ks * km * cf * rhosum / (FaceWidth)))
    )
    Sca = zn * ch * Sc / (kt * kr)
    Mean = Sca - Sigmac
    StdSca = StdScratio * Sca
    StdL = StdLratio * Sigmac
    Std = math.sqrt(StdSca * StdSca + StdL * StdL)
    norm = scipy.stats.norm(Mean, Std)
    rel = 1 - norm.cdf(-Mean / Std)
    return rel


# In[17]:Value assigment shaft
def Shaft_RelVector(data, obj):
    i = -1
    for a in range(0, 29):
        if data.iloc[a]["Name"] == obj.name:
            i = a
    RPM = data.at[i, "RPM_Of_Shaft"]
    Force = data.at[i, "Effective_Force"]
    Torque = data.at[i, "Effective Torque"]
    Length = np.array(obj.Shaft_Length)
    EnduranceLimit = np.array(obj.Shaft_EnduranceLimit)
    UltimateTensileStrength = np.array(obj.Shaft_UltimateTensileStrength)
    ManufactoringProcess = obj.Shaft_ManufactoringProcess
    Diameter = np.array(obj.Shaft_MeanDiameter)
    Temp = data.at[i, "Working_Temperature"]
    Za = np.array(obj.Shaft_EnduranceReliaiblity)
    durn = Shaft_RelLife(
        RPM,
        Force,
        Torque,
        Length,
        EnduranceLimit,
        UltimateTensileStrength,
        ManufactoringProcess,
        Diameter,
        Temp,
        Za,
    )
    YeildStrength = np.array(obj.Shaft_YieldStrength)
    StdYsratio = np.array(obj.Shaft_StdYsratio)
    StdTratio = data.at[i, "StdTratio"]
    StdFratio = data.at[i, "StdFratio"]
    relF = Shaft_Relmoment(
        Force, YeildStrength, Diameter, Length, StdYsratio, StdFratio
    )
    relT = Shaft_Reltorsion(
        Torque, YeildStrength, Diameter, Length, StdYsratio, StdTratio
    )
    res = np.array([relF[0], relT[0], durn[0]])
    return res


# in[17]:value assignment gear
def Gear_RelVector(data_Gear, obj):
    i = -1
    for a in range(0, 28):
        if data_Gear.iloc[a]["Name"] == obj.name:
            i = a
    Load = data_Gear.at[i, "Load"]
    ElasticCoefficient = np.array(obj.Gear_ElasticCoefficient)
    StdScratio = np.array(obj.Gear_StdScratio)
    StdStratio = np.array(obj.Gear_StdStratio)
    StdLratio = data_Gear.at[i, "StdLratio"]
    BSGFactor = np.array(obj.Gear_BSGometryFactor)
    Qualitynumber = data_Gear.at[i, "Qualitynumber"]
    Bhardness = np.array(obj.Gear_BHardness)
    Reliability = np.array(obj.Gear_Reliability)
    Revolution = data_Gear.at[i, "Revolution"]
    LewisFormFactor = np.array(obj.Gear_LewisFormFactor)
    RimThickness = np.array(obj.Gear_RimThickness)
    ToothHight = np.array(obj.Gear_ToothHeight)
    imptdriver = data_Gear.at[i, "ImpactDriver"]
    imptmachine = data_Gear.at[i, "ImpactMachine"]
    PressureAngle = 3.14 * np.array(obj.Gear_PressureAngle) / 180
    NoofTeethp = np.array(obj.Gear_NoofTeeth)
    NoofTeethg = np.array(obj.Gear_NoofTeeth)
    FaceWidth = np.array(obj.Gear_FaceWidth)
    AddendumRadiusg = np.array(obj.Gear_AddendumRadius)
    AddendumRadiusp = np.array(obj.Gear_AddendumRadius)
    DiametralPitchg = np.array(obj.Gear_DiametralPitch)
    DiametralPitchp = np.array(obj.Gear_DiametralPitch)
    Heattreatment = obj.Gear_HeatTreatment
    Material = obj.Gear_Material
    TeethType = obj.Gear_TeethType
    Relc = Gear_RelC(
        NoofTeethg,
        TeethType,
        Material,
        Load,
        ElasticCoefficient,
        StdScratio,
        StdLratio,
        BSGFactor,
        Qualitynumber,
        Bhardness,
        Reliability,
        Revolution,
        LewisFormFactor,
        RimThickness,
        ToothHight,
        imptdriver,
        imptmachine,
        PressureAngle,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
        Heattreatment,
    )
    Relb = Gear_RelB(
        TeethType,
        NoofTeethg,
        Load,
        StdStratio,
        StdLratio,
        BSGFactor,
        Qualitynumber,
        Bhardness,
        Reliability,
        Revolution,
        LewisFormFactor,
        RimThickness,
        ToothHight,
        imptdriver,
        imptmachine,
        PressureAngle,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
        Heattreatment,
        Material,
    )
    life = Gear_life(
        Load,
        PressureAngle,
        NoofTeethg,
        NoofTeethp,
        FaceWidth,
        AddendumRadiusg,
        AddendumRadiusp,
        DiametralPitchg,
        DiametralPitchp,
    )
    res = np.array([Relc[0], Relb[0], life[0]])
    return res


# In[34]:similarity calculator gear
def Gear_similar(data_Gear, onto, obj):
    rel = Gear_RelVector(data_Gear, obj)
    array = []
    for obj1 in onto.Gear.instances():
        if obj != obj1:
            if (
                obj1.Gear_Type == obj.Gear_Type
                and obj1.Gear_Material == obj.Gear_Material
            ):
                a = Gear_RelVector(data_Gear, obj1)
                if (
                    a[0] > (rel[0] - 0.05 * rel[0])
                    and a[0] < (rel[0] + 0.05 * rel[0])
                    and a[1] > (rel[1] - 0.05 * rel[1])
                    and a[1] < (rel[1] + 0.05 * rel[1])
                    and (rel[2] - 0.05 * rel[2]) < a[2]
                    and a[2] < (rel[2] + 0.05 * rel[2])
                ) or (
                    (
                        rel[0] > (a[0] - 0.05 * a[0])
                        and rel[0] < (a[0] + 0.05 * a[0])
                        and rel[1] > (a[1] - 0.05 * a[1])
                        and rel[1] < (a[1] + 0.05 * a[1])
                        and (a[2] - 0.05 * a[2]) < rel[2]
                        and rel[2] < (a[2] + 0.05 * a[2])
                    )
                ):

                    list.append(array, obj1)
    return array


# In[18]:similarity calculator shaft
def Shaft_similar(data_Shaft, onto, obj):
    rel = Shaft_RelVector(data_Shaft, obj)
    array = []
    for obj1 in onto.Shaft.instances():
        if obj != obj1:
            if (
                obj1.Shaft_Type == obj.Shaft_Type
                and obj1.Shaft_Material == obj.Shaft_Material
            ):
                a = Shaft_RelVector(data_Shaft, obj1)
                if (
                    a[0] > rel[0] - 0.05 * rel[0]
                    and a[0] < rel[0] + 0.05 * rel[0]
                    and a[1] > rel[1] - 0.05 * rel[1]
                    and a[0] < rel[1] + 0.05 * rel[1]
                    and a[2] > rel[2] - 0.05 * rel[2]
                    and a[2] < rel[2] + 0.05 * rel[2]
                ):
                    list.append(array, obj1)
    return array


# In[20]: main function to take input from user and call appropriate functions
print(
    "Available Machines ID is given by MachID where ID can take value from 1 to 28. Eg Mach3 "
)
print("Please enter the name of the machine from above")
Machine_ID = input()
print("Please enter the type of Component (Shaft\ Gear)")
Type = input()
TTFS = "TTFShaft.xlsx"
TTFG = "TTFGear.xlsx"
df = pd.DataFrame()
A = np.array([])
Y = np.array([])
Mach = []
for obj in onto.Machine.instances():
    if obj.name == Machine_ID:
        Mach = obj.Has_Components
if Mach == []:
    print("No Machine found with that name")
    sys.exit(1)
for obj1 in Mach:
    if Type == "Shaft":
        for obj2 in onto.Shaft.instances():
            if obj1 == obj2:
                df = pd.read_excel(TTFS, obj2.name)
                sim = Shaft_similar(data_Shaft, onto, obj2)
                if sim == []:
                    print("No similar Shaft present on other machine in Network")
                    sys.exit(1)
                else:
                    print("Similar Shaft are present on given Machine.")
                for shaft in sim:
                    obje = shaft.ComponentOf
                    print(obje[0].name)
                    df1 = pd.read_excel(TTFS, shaft.name)
                    df = df.append(df1, ignore_index=True)

    elif Type == "Gear":
        for obj2 in onto.Gear.instances():
            if obj1 == obj2:
                df = pd.read_excel(TTFG, obj2.name)
                Sim = Gear_similar(data_Gear, onto, obj2)
                if Sim == []:
                    print("No similar Gear present on other machine in Network ")
                    sys.exit(1)
                else:
                    print("Similar Gear are present on given Machine.")
                for gear in Sim:
                    objec = gear.ComponentOf
                    print(objec[0].name)
                    df1 = pd.read_excel(TTFG, gear.name)
                    df = df.append(df1, ignore_index=True)
    else:
        print("Given type is not added in the system. We are working on it")
count = len(df)
A = np.array(df["TimeToFailureA"])
for i in range(count):
    Y = np.append(
        Y,
        reliability(
            df.at[i, "TimeToFailureA"], df.at[i, "Beta"], df.at[i, "TimeToFailureT"]
        ),
    )
print(
    "Please enter the reliability of the component for which you want to compute working Hours"
)
Tn = float(input())
a = linearreg(A, Y, Tn)
print("The selected component will have", Tn, "reliability at", a, "working Hours")
