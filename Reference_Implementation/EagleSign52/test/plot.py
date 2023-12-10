# Implementors: EagleSign Team
# This implementation is highly inspired from Dilithium and
# Falcon Signatures' implementations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

LEVEL = "52"


if __name__=="__main__":
    # Fg
    data = pd.read_csv("distributionFg_EagleSign%s.csv"%(LEVEL), sep="\t")
    top = max(list(data["count"])[1:])
    right =list(data.coef)[-1]

    fig = plt.figure()
    ax = fig.add_subplot()
    fig.subplots_adjust(top=0.85)

    fig.suptitle('EagleSign %s: Probability Distribution of Fg'%(LEVEL), fontsize=14, fontweight='bold')
    ax.set_title('$\mathbb{P}(X = k), k \in [%s, %s]$'%(list(data.coef)[1],list(data.coef)[-1]))
    ax.text(right - 3.5*right/7, top - 0.03*top/0.175, '$\overline{X}= %5f$\n$\sigma= %5f$'%(list(data.coef)[0], np.sqrt(list(data["count"])[0])), style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})

    ax.set_xlabel('$k \in [%s, %s]$'%(list(data.coef)[1],list(data.coef)[-1]))
    ax.set_ylabel('$\mathbb{P}(X = k)$')
    ax.plot(list(data.coef)[1:], list(data["count"])[1:])
    plt.savefig("EagleSign%sFg.png"%(LEVEL))


    # DF
    data = pd.read_csv("distributionDF_EagleSign%s.csv"%(LEVEL), sep="\t")
    top = max(list(data["count"])[1:])
    right =list(data.coef)[-1]

    fig = plt.figure()
    ax = fig.add_subplot()
    fig.subplots_adjust(top=0.85)

    fig.suptitle('EagleSign %s: Probability Distribution of DF'%(LEVEL), fontsize=14, fontweight='bold')
    ax.set_title('$\mathbb{P}(X = k), k \in [%s, %s]$'%(list(data.coef)[1],list(data.coef)[-1]))
    ax.text(right - 3.5*right/7, top - 0.03*top/0.175, '$\overline{X}= %5f$\n$\sigma= %5f$'%(list(data.coef)[0], np.sqrt(list(data["count"])[0])), style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})

    ax.set_xlabel('$k \in [%s, %s]$'%(list(data.coef)[1],list(data.coef)[-1]))
    ax.set_ylabel('$\mathbb{P}(X = k)$')
    ax.plot(list(data.coef)[1:], list(data["count"])[1:])
    plt.savefig("EagleSign%sDF.png"%(LEVEL))
