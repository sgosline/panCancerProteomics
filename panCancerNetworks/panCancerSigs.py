# -*- coding: utf-8 -*-
"""
Goal is to create a single set of pancancer signatures

@author: gosl241
"""

import argparse
import sys

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pandas as pd
import pickle
import random

class kvdictAppendAction(argparse.Action):
    """
    argparse action to split an argument into KEY=VALUE form
    on the first = and append to a dictionary.
    """
    def __call__(self, parser, args, values, option_string=None):
        assert(len(values) == 1)
        values = values[0].split(',')
        for val in values:
            try:
                (k, v) = val.split("=", 2)
            except ValueError as ex:
                raise argparse.ArgumentError(self, \
                                             "could not parse argument \"{values[0]}\" as k=v format")
            d = getattr(args, self.dest) or {}
            d[k] = v
            setattr(args, self.dest, d)


#Parser information for command line
parser = argparse.ArgumentParser(description="""Get data from the proteomic \
                                 data commons and build community networks""")
parser.add_argument('--enrich', dest='doEnrich', action='store_true',\
                    default=False, help='Flag to do GO enrichment')
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')
parser.add_argument('--fromFile',dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')
parser.add_argument('--quantile',dest='qt',default=0.01,type=float,\
                    help='Threshold to select top proteins from each patient')
parser.add_argument('--sample', dest='sample', default=False, action='store_true',\
                    help='Use this flag if you want to sample 5 patients from each disease to test')
parser.add_argument('--graph', dest='graph', default='../../data/igraphPPI.pkl',\
                help='Path to pickled igraph interactome')


def loadCancerData(qt):
    norms = prot.normals_from_manifest('../cancerProtFeatures/data/PDC_biospecimen_manifest_07182020_151323.csv')

#    bcData = prot.parsePDCfile('../cancerProtFeatures/data/TCGA_Breast_BI_Proteome.itraq.tsv')
    bcData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')
    hccData = prot.parsePDCfile('../cancerProtFeatures/data/Zhou_Liver_Cancer_Proteome.tmt11.tsv')
    hnsccData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC3_Head_and_Neck_Carcinoma_Proteome.tmt11.tsv')
    ovcaData = prot.parsePDCfile('../cancerProtFeatures/data/TCGA_Ovarian_PNNL_Proteome.itraq.tsv')

    normPats = {'brca': set([a for a in bcData['Patient'] if a in norms['Breast Invasive Carcinoma']]),\
                'coad': set([a for a in colData['Patient'] if a in norms['Colon Adenocarcinoma']]),\
                'luad': set([a for a in lungData['Patient'] if a in norms['Lung Adenocarcinoma']]),\
                'gbm': set([a for a in gbmData['Patient'] if a in norms['Other']]),
                'hcc': set([a for a in hccData['Patient'] if a in ['P'+n for n in norms['Hepatocellular Carcinoma ']]]),\
                'ovca': set([a for a in ovcaData['Patient'] if a in norms['Ovarian Serious Cystadenocarcinoma']]),\
                'hnscc': set([a for a in hnsccData['Patient'] if a in norms['Head and Neck Squamous Cell Carcinoma']])}

    for key, np in normPats.items():
        print(len(np), 'normals for', key)

    namemapper = None #hyp.mapHGNCtoNetwork()

    #here we get the top most distinguished from normals
    patDiffs = {'hnscc': prot.getTumorNorm(hnsccData, normPats['hnscc'], namemapper, quantThresh=qt),
               # 'ovca': prot.getTumorNorm(ovcaData, normPats['ovca'], namemapper, quantThresh=qt),
                'hcc': prot.getTumorNorm(hccData, normPats['hcc'], namemapper, quantThresh=qt),
                'brca': prot.getTumorNorm(bcData, normPats['brca'], namemapper, quantThresh=qt),
                'luad': prot.getTumorNorm(lungData, normPats['luad'], namemapper, quantThresh=qt),
                'coad': prot.getTumorNorm(colData, normPats['coad'], namemapper, quantThresh=qt),
                'gbm': prot.getTumorNorm(gbmData, normPats['gbm'], namemapper, quantThresh=qt)}
    pc = dict()
    for key,val in patDiffs.items():
        pc.update(val)
    print("PanCan dictionary has",len(pc),'patients')
    patDiffs={'panCan':pc}
    return patDiffs

def build_hyphae_from_data(qt, g, sample=False):
    """ Temp function to load data from local directory"""
    ##this is the framework for the PDC data parser.
    #now we want to build network communities for each
    hyphae = dict()
    patDiffs = loadCancerData(qt)
    beta = 0.5
    for key, vals in patDiffs.items():
        if sample:
            new_vals = {}
            for v in random.sample(list(vals), 300):
                new_vals[v] = vals[v]
            vals = new_vals
            print(len(vals))
        this_hyp = hyphalNetwork(vals, g.copy(),beta=beta, g=3, do_forest=False, noComms=False)
        hyphae[key+str(qt)] = this_hyp
        this_hyp._to_file(key+str(qt)+'_hypha.pkl')
    return hyphae

def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae


def main():
    args = parser.parse_args()
    gfile = args.graph
    g = pickle.load(open(gfile, 'rb'))#hyp.make_graph_from_dict(gfile)
    if args.fromFile is None:
        hyphae = build_hyphae_from_data(args.qt, g, args.sample)
    else:
        hyphae = loadFromFile(args.fromFile)

    #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = hyStats.compute_all_distances(hyphae)
        res.to_csv('panCancerDistances.csv')
        nmi = hyStats.compute_all_nmi(hyphae, g)
        nmi.to_csv('panCancerNMI.csv')

    for key, this_hyp in hyphae.items():
        this_hyp.node_stats().to_csv(key+'_nodelist.csv')
        if args.doEnrich:
            if len(this_hyp.community_enrichment) == 0:
                com_e = hyEnrich.go_enrich_communities(this_hyp)
                this_hyp.assign_enrichment(com_e, type='community')
                this_hyp._to_file(key+'_hypha.pkl')
                com_e.to_csv(key+'enrichedCommunityGOterms.csv')
            ##next: compare enrichment between patients mapped to communities
        this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')
        this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')

if __name__ == '__main__':
    main()

if True:
    g = pickle.load(open('../../data/igraphPPI.pkl', 'rb'))
    hyphae = build_hyphae_from_data(0.01, g, True)
