import json
from os.path import exists
from boolean3_addon import attr_cy
from sbie_aging.results import table_s1
from os.path import dirname

analysisresult = dirname(table_s1.__file__)+'/c/attractors.json'

modeltext = '''
RAS = Random
Rheb = Random
IRS = Random
SIRT1 = Random
ERK = Random
G6pase = Random
p38MAPK = Random
ULK1 = Random
FOUR_EBP1 = Random
FOXO = Random
IGF1R = Random
ROS = Random
S6K = Random
PGC1a = Random
ATG13 = Random
PDK = Random
PP2A = Random
MKK = Random
PTEN = Random
AKT = Random
LC = Random
MnSOD = Random
AXP = Random
PI3K = Random
p16 = Random
Bcl2 = Random
CDK = Random
DNAdamage = Random
LowNuEx = Random
AMPK = Random
IGF1 = Random
BIM = Random
BAX = Random
ATG5 = Random
TSC = Random
pRB1 = Random
caspase3 = Random
p53 = Random
mTOR = Random
E2F = Random
p21 = Random
NAD = Random
p27 = Random
ATP = Random
SGK = Random
Glycolysis = Random
BMI1 = Random
TCAcyc = Random
GLUT4 = Random
MDM2 = Random
Glucose = Random
AKT *= (PDK) and not (PP2A)
AMPK *= (AXP or p53)
ATG13 *=  not (mTOR)
ATG5 *= (FOXO)
ATP *= (Glycolysis or TCAcyc)
AXP *= (LowNuEx) and not (ATP)
BAX *= (p53 or BIM) and not (AKT or Bcl2)
BIM *= (FOXO) and not (AKT or ERK)
BMI1 *= (AKT)
Bcl2 *= (AKT) and not (p53)
CDK *=  not (p16 or p21)
DNAdamage *= (ROS)
E2F *=  not (pRB1)
ERK *= (RAS) and not (PP2A)
FOUR_EBP1 *=  not (mTOR)
FOXO *= (AMPK or PP2A or SIRT1) and not (AKT or ERK or MDM2 or SGK)
G6pase *= (FOXO)
GLUT4 *= (AMPK or PGC1a)
Glucose *= (GLUT4)
Glycolysis *= (G6pase or Glucose)
IGF1R *= (IGF1)
IRS *= (IGF1R) and not (S6K)
LC *= (ATG13 or ATG5 or ULK1)
MDM2 *= (AKT) and not (S6K)
MKK *= (DNAdamage)
MnSOD *= (FOXO)
NAD *= (AMPK)
PDK *= (PI3K)
PGC1a *= (AMPK or mTOR or SIRT1)
PI3K *= (IRS) and not (PTEN)
PP2A *=  not (mTOR)
PTEN *= (p53)
RAS *= (IGF1R)
ROS *= (TCAcyc) and not (MnSOD)
Rheb *=  not (TSC)
S6K *= (mTOR)
SGK *= (p53 or PDK)
SIRT1 *= (FOXO or NAD)
TCAcyc *= (Glycolysis)
TSC *= (AMPK or SIRT1) and not (AKT or ERK)
ULK1 *=  not (AMPK or mTOR)
caspase3 *= (BAX)
mTOR *= (Rheb) and not (AMPK)
p16 *= (p38MAPK) and not (BMI1)
p21 *= (FOXO or p53) and not (AKT)
p27 *= (FOXO) and not (AKT)
p38MAPK *= (MKK) and not (PP2A)
p53 *= (AMPK or DNAdamage or p38MAPK) and not (MDM2 or SIRT1)
pRB1 *=  not (CDK)
'''

attr_cy.build(modeltext)

import pyximport; pyximport.install()

res = attr_cy.run(samples=5000, steps=1000, 
    debug=False, on_states=[], progress=True)

json.dump(res, open(analysisresult, 'w'), indent=4)


