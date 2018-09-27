import itertools
import numpy as np

def generate_twonode_topologies() :
    options = [-1, 0, 1]

    topology_list, logic_list = [], []

    for [AA, AB, BA, BB] in itertools.product(options, repeat=4) :
        if BA == 0 :
            continue

        [FA, FB] = [0, 0]
        if AA >= 0 and BA >= 0 :
            FA = -1
        if not (AB == 0 and BB == 0) :
            if (AB >= 0 and BB >= 0) :
                FB = -1
            elif (AB <= 0 and BB <= 0) :
                FB = 1

        top = [AA, AB, FA, BA, BB, FB]

        logics = ['AorBor']
        if (AA > 0) or (BA > 0) or (AA != 0 and AA == BA) :
            logics.append('AandBor')
        if BB != 0 and BB == AB :
            logics.append('AorBand')
        if ((AA > 0) or (BA > 0) or (AA != 0 and AA == BA)) and BB != 0 and BB == AB :
            logics.append('AandBand')

        #for logic in logics :
        #    print top, logic

        topology_list.append(top)
        logic_list.append(logics)

    return topology_list, logic_list

def generate_threenode_topologies() :
    options = [-1, 0, 1]

    topology_list, logic_list = [], []

    for [AA, AB, AC, BA, BB, BC, CA, CB, CC] in itertools.product(options, repeat=9) :
        # if signal can't get to output node
        if (AC == 0) and (AB == 0 or BC == 0) :
            continue
        # if B doesn't do anything
        if BA == 0 and BC == 0 :
            continue

        [FA, FB, FC] = [0, 0, 0]
        if AA >= 0 and BA >= 0 and CA >= 0 :
            FA = -1
        if not (AB == 0 and BB == 0 and CB == 0) :
            if (AB >= 0 and BB >= 0 and CB >= 0) :
                FB = -1
            elif (AB <= 0 and BB <= 0 and CB <= 0) :
                FB = 1
        if not (AC == 0 and BC == 0 and CC == 0) :
            if (AC >= 0 and BC >= 0 and CC >= 0) :
                FC = -1
            elif (AC <= 0 and BC <= 0 and CC <= 0) :
                FC = 1

        top = [AA, AB, AC, FA, BA, BB, BC, FB, CA, CB, CC, FC]

        logics = ['AorBorCor']
        Aand, Band, Cand = False, False, False
        if (AA > 0 or BA > 0 or CA > 0) or (sum([int(x < 0) for x in [AA, BA, CA]]) >= 2) :
            Aand = True
        if (sum([int(x > 0) for x in [BB, AB, CB]]) >= 2) or (sum([int(x < 0) for x in [BB, AB, CB]]) >= 2) :
            Band = True
        if (sum([int(x > 0) for x in [CC, AC, BC]]) >= 2) or (sum([int(x < 0) for x in [CC, AC, BC]]) >= 2) :
            Cand = True

        if Aand :
            logics.append('AandBorCor')
        if Band :
            logics.append('AorBandCor')
        if Cand :
            logics.append('AorBorCand')
        if Aand and Band :        
            logics.append('AandBandCor')
        if Aand and Cand :        
            logics.append('AandBorCand')
        if Band and Cand :        
            logics.append('AorBandCand')
        if Aand and Band and Cand :        
            logics.append('AandBandCand')

        #for logic in logics :
        #    print top, logic

        topology_list.append(top)
        logic_list.append(logics)

    return topology_list, logic_list


def sample_parameters(topology, num_samples=1000) :

    # inputs:
    # topology: an array indicating the presence and sign of a topology link
    # for example, [1, -1, -1, 0]
    # num_samples: number of parameter samples to generate

    # return:
    # a list of lists of tuples
    # each tuple pair is (kcat, Km) for one topology link
    # inner list: samples for each link
    # outer list: topology links

    # sign of each topology link is encoded in the sign of kcat

    kcat = { 'min' : 0.1, 'max' : 10, 'samples' : [] }
    Km = { 'min' : 0.001, 'max' : 1000, 'samples' : [] }

    kcat['values'] = np.logspace(np.log10(kcat['min']), np.log10(kcat['max']), num_samples)
    Km['values'] = np.logspace(np.log10(Km['min']), np.log10(Km['max']), num_samples)

    for link in topology :
        if link == 0 :
            kcat_samples = [0]*num_samples
            Km_samples = [0]*num_samples
        else :
            kcat_samples = [x*link for x in np.random.choice(kcat['values'], size=num_samples, replace=False)]
            Km_samples = list(np.random.choice(Km['values'], size=num_samples, replace=False))

        kcat['samples'].append(kcat_samples)
        Km['samples'].append(Km_samples)

    return [zip(x,y) for (x,y) in zip(kcat['samples'], Km['samples'])]


if __name__ == '__main__' :

    top2, logic2 = generate_twonode_topologies()
    top3, logic3 = generate_threenode_topologies()

    # example: sample parameters for first topology in the 2-node set
    # generate 5 parameter samples
    print sample_parameters(top2[0], num_samples=5)
    