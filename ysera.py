import pandas as pd
import numpy as np
import math
import os
from os import listdir
from os.path import isfile, join
from sklearn.metrics.pairwise import euclidean_distances
import time
import gzip
import json
PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))

HB_DEFAULT = 3.1
SB_DEFAULT = 4.0
DB_DEFAULT = 2.2
VDW_DEFAULT = 5.5
PS_DEFAULT = 7.2
AAAN_BEG_DEFAULT = 2.0
AAAN_END_DEFAULT = 3.0
AASPI_DEFAULT = 5.3
AACTN_BEG_DEFAULT = 3.4
AACTN_END_DEFAULT = 4.0

def ysera(filename, params):
    
    if(not("hb" in params)):
        params['hb'] = HB_DEFAULT
    if(not("sb" in params)):
        params['sb'] = SB_DEFAULT
    if(not("db" in params)):
        params['db'] = DB_DEFAULT
    if(not("vdw" in params)):
        params['vdw'] = VDW_DEFAULT
    if(not("ps" in params)):
        params['ps'] = PS_DEFAULT
    if(not("aaan_beg" in params)):
        params['aaan_beg'] = AAAN_BEG_DEFAULT
    if(not("aaan_end" in params)):
        params['aaan_end'] = AAAN_END_DEFAULT
    if(not("aaspi" in params)):
        params['aaspi'] = AASPI_DEFAULT
    if(not("aactn_beg" in params)):
        params['aactn_beg'] = AACTN_BEG_DEFAULT
    if(not("aactn_end" in params)):
        params['aactn_end'] = AACTN_END_DEFAULT

    print(params)
    
    res = myfunction(filename, params)
    return res 

def myfunction(filename, params):
    string1 = ""
    string2 = ""
    
    path = PROJECT_HOME + '/temp/' + filename
    pathoutput = PROJECT_HOME + '/output/' + filename

    AROMTRP = ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
    AROMPHE = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    AROMTYR = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    
    datapdb = pd.DataFrame(
        columns=['Type', 'Atom ID', 'Atom Type', 'aa', 'Chain ID', 'X', 'Y', 'Z', 'occupancy', 'temperature factor',
                    'element symbol'])
    i = 0
    
    print("Starting to load data")
    start_time = time.time()
    Amin = ''
    Aromaticpos = []
    AromaticArray = {}
    AromaticPoints = []
    AromaticNormals = {}
    Invalids = []
    
    file = open(path, 'r')
    
    print(file)
    
    for line in file:
        if ("ENDMDL" in line):
            break
        if (line[16:20].strip() in "HOH"):
            continue
        if ("ATOM" in line[:6] or "HETATM" in line[:6]):

            if (line[17:20].strip() in ['TYR', 'PHE', 'TRP']):
                # print("{}{}".format(Amin,line[17:20].strip()))
                if (Amin == ''):
                    Amin = line[20:27].strip()
                elif (Amin != line[20:27].strip()):
                    AromaticArray[Amin] = np.asarray(Aromaticpos)
                    if (len(AromaticPoints) == 3):
                        veca = np.subtract(AromaticPoints[1], AromaticPoints[0])
                        vecb = np.subtract(AromaticPoints[2], AromaticPoints[0])
                        AromaticNormals[Amin] = np.cross(veca, vecb)
                    else:
                        print("Incomplete Aromatic Structure at {}".format(Amin))
                        Invalids.append(Amin)
                    Amin = line[20:27].strip()
                    Aromaticpos = []
                    AromaticPoints = []
                elif ((line[17:20].strip() == 'TYR' and line[11:16].strip() in AROMTYR) or (
                        line[17:20].strip() == 'PHE' and line[11:16].strip() in AROMPHE) or (
                                line[17:20].strip() == 'TRP' and line[11:16].strip() in AROMTRP)):
                    if (Aromaticpos == []):
                        Aromaticpos = [float(line[27:38].strip()), float(line[38:46].strip()),
                                        float(line[46:54].strip())]
                    else:
                        Aromaticpos = [(x + y) / 2 for x, y in zip(Aromaticpos, [float(line[27:38].strip()),
                                                                                    float(line[38:46].strip()),
                                                                                    float(line[46:54].strip())])]
                    if (len(AromaticPoints) < 3):
                        AromaticPoints.append(
                            [float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
            datapdb.loc[i] = [line[:6].strip(), line[6:11].strip(), line[11:16].strip(), line[17:20].strip(),
                                line[20:27].strip(), line[27:38].strip(), line[38:46].strip(), line[46:54].strip(),
                                line[54:60].strip(), line[60:66].strip(), line[66:80].strip()]
            # datapdb.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
            i += 1
    AromaticArray[Amin] = np.asarray(Aromaticpos)
    if (len(AromaticPoints) == 3):
        veca = np.subtract(AromaticPoints[1], AromaticPoints[0])
        vecb = np.subtract(AromaticPoints[2], AromaticPoints[0])
        AromaticNormals[Amin] = np.cross(veca, vecb)
    elif (0 < len(AromaticPoints) < 3):
        print("Incomplete Aromatic Structure at {}".format(Amin))
        Invalids.append(Amin)
    print("DataLoaded")
    print("---%s seconds ---" % (time.time() - start_time))
    if 'Atom ID' in datapdb.columns:
        datapdb = datapdb.drop(['Atom ID', 'occupancy', 'temperature factor', 'element symbol'], axis=1)
    datapdb["X"] = np.float32(datapdb["X"], downcast='signed')
    datapdb["Y"] = np.float32(datapdb["Y"], downcast='signed')
    datapdb["Z"] = np.float32(datapdb["Z"], downcast='signed')

    dist = euclidean_distances(np.float32(datapdb[["X", "Y", "Z"]].to_numpy()),
                                np.float32(datapdb[["X", "Y", "Z"]].to_numpy()))
    Dist = pd.DataFrame(data=dist, index=None)

    New = pd.concat([datapdb, Dist], axis=1, sort=False)

    hb = 0
    lighbacep = ['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1']  # 3
    lighbdono = ['HE1', 'H', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HD1', 'HE2', 'HG', 'HG1', 'HG21', 'HG22',
                    'HG23', 'HH', 'HD21', 'HD22', 'HE21', 'HE22']  # 18

    sb = 0
    ligsb1 = ['OD2', 'OE2', 'OH']
    ligsb2 = ['NZ', 'NH2', 'NH1']
    aasbneg = ['ASP', 'GLU']
    aasbpos = ['LYS', 'ARG']

    db = 0
    ligdb = ['SG']

    vdw = 0
    ligvdw = ['CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE']
    aavdw = ['VAL', 'TRE', 'MET', 'LEU', 'ILE']

    lpi = 0
    tshaped = 0
    inter = 0
    paralel = 0
    aapi = ['TYR', 'PHE', 'TRP']

    ctn = 0
    ligctn = ['MG', 'CU', 'K', 'FE2', 'FE', 'NI', 'NA', 'MO1', 'MO3', 'MO4', 'MO5', 'MO6', 'MO7', 'MO8', 'MO9',
                'NZ', 'NH2', 'NH1']
    ligctn2 = ['CG', 'CE2', 'CG']
    aactn = ['TYR', 'PHE', 'TRP']

    an = 0
    ligan1 = ['CL', 'BR', 'I', 'OD2', 'OE2', 'OH']
    ligan2 = ['CG', 'CE2', 'CG']
    aaan = ['TYR', 'PHE', 'TRP']

    spi = 0
    ligspi1 = ['SG']
    ligspi2 = ['CG', 'CE2', 'CG']
    aaspi = ['TYR', 'PHE', 'TRP']

    f = open('output/' + filename + ".txt", "w+")
    print("Starting to process data")
    start_time = time.time()
    Exclusions = []
    for i in range(len(New)):
        for j in range(i + 1, len(New)):
            distance = New[j].iloc[i]
            if (distance > 8 or distance == 0):
                continue
            atom1 = New['Atom Type'].iloc[i]
            atom2 = New['Atom Type'].iloc[j]
            aa1 = New['aa'].iloc[i]
            aa2 = New['aa'].iloc[j]
            chaincode1 = New['Chain ID'].iloc[i]
            chaincode2 = New['Chain ID'].iloc[j]
            distance = New[j].iloc[i]
            # print(chaincode1,chaincode2)
            # HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---
            if atom1 in lighbacep[:] and atom2 in lighbdono[:] and 0.0 < distance < params['hb']:
                hb += 1
                string2 = (
                    'Hydrogen_Bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            elif atom1 in lighbdono[:] and atom2 in lighbacep[:] and 0.0 < distance < params['hb']:
                hb += 1
                string2 = string2 + (
                    'Hydrogen_Bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            # SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----
            if aa1 in aasbpos[:] and atom1 in ligsb2[:] and aa2 in aasbneg[:] and atom2 in ligsb1[
                                                                                            :] and 0.0 < distance < params['sb']:
                sb += 1
                string2 = string2 + (
                    'Salt_Bridge' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            elif aa1 in aasbneg[:] and atom1 in ligsb1[:] and aa2 in aasbpos[:] and atom2 in ligsb2[
                                                                                                :] and 0.0 < distance < params['sb']:
                sb += 1
                string2 = string2 + (
                    'Salt_Bridge' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            # DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----
            if atom1 in ligdb[:] and atom2 in ligdb[:] and 0.0 < distance < 2.2:
                db += 1
                string2 = string2 + (
                    'Dissulfide_bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            # VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------
            if aa1 in aavdw[:] and atom1 in ligvdw[:] and aa2 in aavdw[:] and atom2 in ligvdw[
                                                                                        :] and 0.0 < distance < params['vdw']:
                vdw += 1
                print('{}'.format(distance))
                string2 = string2 + (
                    'Van_der_Waals' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                        distance) + '\n')
            # PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------
            if aa1 in aapi[:] and aa2 in aapi[:]:
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = AromaticArray[chaincode1]
                    coordinates2 = AromaticArray[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (aromaticdistance < params['aaspi']):
                        string2 = string2 + (
                            'Pi_stacking  ' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        NormalVector1 = AromaticNormals[chaincode1] / np.linalg.norm(AromaticNormals[chaincode1])
                        NormalVector2 = AromaticNormals[chaincode2] / np.linalg.norm(AromaticNormals[chaincode2])
                        Angle = np.arccos(np.clip(np.dot(NormalVector1, NormalVector2), -1.0, 1.0))
                        if (Angle > 50):
                            tshaped += 1
                        elif (30 < Angle < 50):
                            inter += 1
                        elif (Angle < 30):
                            paralel += 1
                        # print(Angle,type(AromaticNormals[chaincode1]),NormalVector1,NormalVector2,np.dot(NormalVector1, NormalVector2))
                        lpi += 1
                        Exclusions.append([chaincode1, chaincode2])
            # CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------
            if aa1 in aactn[:] and atom2 in ligctn[:]:
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = AromaticArray[chaincode1]
                    coordinates2 = np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]])
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aactn_beg'] < aromaticdistance < params['aactn_end']):
                        ctn += 1
                        string2 = string2 + (
                            'Cation_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            elif (atom1 in ligctn[:] and aa2 in aactn[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = AromaticArray[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aactn_beg'] < aromaticdistance < params['aactn_end']):
                        ctn += 1
                        string2 = string2 + (
                            'Cation_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            # SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------
            if (aa1 in aaspi[:] and atom2 in ligspi1[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = AromaticArray[chaincode1]
                    coordinates2 = np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]])
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (0 < aromaticdistance < params['aaspi']):
                        spi += 1
                        string2 = string2 + (
                            'Sulfur_Aryl  ' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            elif (atom1 in ligspi1[:] and aa2 in aaspi[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = AromaticArray[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (0 < aromaticdistance < params['aaspi']):
                        spi += 1
                        string2 = string2 + (
                            'Sulfur_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            # AN-------------AN-------------AN-------------AN-------------AN-------------AN-------------
            if aa1 in aaan[:] and atom2 in ligan1[:]:
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = AromaticArray[chaincode1]
                    coordinates2 = np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]])
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aaan_beg'] < aromaticdistance < params['aaan_end']):
                        an += 1
                        string2 = string2 + (
                            'Anion_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            elif (atom1 in ligan1[:] and aa2 in aaan[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,
                                                                                    chaincode2] not in Exclusions and [
                    chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = AromaticArray[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aaan_beg']< aromaticdistance < params['aaan_end']):
                        an += 1
                        string2 = string2 + (
                            'Anion_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                                aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])


    f.write(string2)
    f.close()
    string1 = { "filename" : filename, 
        "hb" : hb,
        "sb" : sb,
        "db" : db,
        "lpi" : lpi,
        "tshaped": tshaped,
        "inter" : inter,
        "paralel": paralel,
        "vdw": vdw,
        "ctn": ctn,
        "an": an,
        "spi": spi
    }
    return string1
 