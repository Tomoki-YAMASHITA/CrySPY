from logging import getLogger

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

from ..IO import read_input as rin
from ..IO import pkl_data


logger = getLogger('cryspy')
"""
c_struc : current gen structure {ID:Energy..}   eV/atomとなっているので注意
vc_nat : {ID:nat..}  vcc_nat                ea_dataの読み込みのようにする
end_point : [end_point0,end_point1...] rin.end_point

rsltのエネルギーが1原子あたりかを確認

"""


def calc_ef(energy, ratio, end_point):
    '''
    energy: eV/atom
    
    note:
    nat = [4, 3, 5]
    12 atoms: energy * 12 - 4 * end_point[0] - 3 * end_point[1] - 5 * end_point[2]
    1 atom: energy - 4/12 * end_point[0] - 3/12 * end_point[1] - 5/12 * end_point[2]
    --> energy - ratio[0] * end_point[0] - ratio[1] * end_point[1] - ratio[2] * end_point[2]
    '''
    # ---------- np.nan
    if np.isnan(energy):
        return np.nan
    # ---------- check
    if len(ratio) != len(end_point):
        logger.error('len(ratio) != len(end_point)')
        raise SystemExit(1)
    # ---------- calc formation energy
    ef = energy
    for x, y in zip(ratio, end_point):
        ef -= x * y
    # ---------- return
    return ef


def draw_convex_hull(vpoints, x, ef, c_struc_num):
    # ---------- draw setting
    # ---------- figure
    plt.rcParams['figure.figsize'] =[8, 6]
    plt.rcParams["figure.dpi"] = 120
    plt.rcParams['figure.facecolor'] = 'white'

    # ---------- axes
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.linewidth'] = 1.5

    # ---------- ticks
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['xtick.major.size'] = 8.0
    plt.rcParams['ytick.major.size'] = 8.0

    # ---------- lines
    plt.rcParams['lines.linewidth'] = 2.5
    plt.rcParams['lines.markersize'] = 12

    # ---------- grid
    plt.rcParams['grid.linestyle'] = ':'

    # ---------- font
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.labelsize'] = 26
    plt.rcParams['legend.fontsize'] = 26
    #plt.rcParams['pdf.fonttype'] = 42    # embed fonts in PDF using type42 (True type)

    if(max(ef) > 0):
        max_ef = 0.05

    fig, ax = plt.subplots(1, 1)
    ax.axhline(y=0, xmin=0, xmax=1, color='black', linestyle='--')
    ax.plot(vpoints[1:, 0], vpoints[1:, 1], color='C1')
    ax.set_ylim(min(ef)*1.8,max_ef)
    
    # plot red:current generation / gray:other generation
    if c_struc_num == len(ef):
        ax.scatter(x,ef,c='orange')
    else:
        for i in range(len(ef) - c_struc_num):
            ax.scatter(x[i],ef[i],c='gray')

        for i in range(c_struc_num):
            ax.scatter(x[-i-1],ef[-i-1],c='orange')
            logger.info('ax.scatter_x:'+str(x[-i-1]))
            logger.info('ax.scatter_y:'+str(ef[-i-1]))


def write_asc_hdist(hdist):
    (gen, id_queueing,id_running) = pkl_data.load_ea_id()
    with open('./data/asc_hdist','a') as asc_h:
        asc_h.write(f"----- Gen{gen} -----\n")
        asc_h.write('ID'+'\t'+'hull distance'+'\n')
        asc_h.writelines("\n".join(str(i)+"\t"+str(j) for i,j in hdist.items()))
        asc_h.write('\n\n')
        

    
    


def cal_hull_distance2d(xdata,ydata,equations):
    hullE = []
    for i,equation in enumerate(equations):
        #logger.info(str(xdata)+' '+str(ydata))
        hullE.append(abs(ydata-(-1*equation[0]*xdata-1*equation[2])/equation[1]))
    return min(hullE)
    



def calc_convex_hull(ef,c_struc_num):
    """
    ef: 全ての構造のefの辞書、rslt_data.pklから持ってきた {ID: Energy...}
    x: EAvc_data.pklから持ってきた割合
    c_struc_num: number of current generation sturc 
    エリート構造選択時にも呼び出す
    """

    # ---------- calc convex hull and hull distance
    vc_nat,x = pkl_data.load_eavc_data()
    rate_x = []
    rslt_ef = []
    logger.info('x is'+str(x))
    logger.info('ef is'+str(ef))

    # ------ convert dic to list
    """
    エリート構造選択時に呼び出した際efは現在の構造までしかないが
    xは次の世代のデータも入っているのでefの数だけ取り出す
    """
    for i in range(len(ef)):
        rate_x.append(x[i])
        rslt_ef.append(ef[i])

    # ------ add end_point value
    rate_x.insert(0,1)
    rate_x.insert(0,0)
    rslt_ef.insert(0,0)
    rslt_ef.insert(0,0)

    #check x and ef
    # logger.info('x is'+str(rate_x))
    # logger.info('ef is'+str(rslt_ef))

    # ------ get for plot value
    plot_x = []
    plot_ef = []
    for i in range(len(rslt_ef)):
        if rslt_ef[i] <= 0:
            plot_x.append(rate_x[i])
            plot_ef.append(rslt_ef[i])

    # ------ calc convex hull
    #efがマイナスのものだけをプロットするようにしているので場合によってはマイナスのefが少なすぎてプロットできないかもしれない
    #改良が必要
    try:
        points=list(zip(plot_x,plot_ef))
        hull = ConvexHull(points)
    except Exception as e:
        points=list(zip(rate_x,rslt_ef))
        hull = ConvexHull(points)

    #logger.info('hull_points:'+str(hull.points))
    vpoints = hull.points[hull.vertices]
    vpoints = np.vstack((vpoints, vpoints[0]))

    #logger.info('equations:'+str(hull.equations))

    # ------ draw convex hull
    draw_convex_hull(vpoints,rate_x,rslt_ef,c_struc_num)
    logger.info(str(len(rate_x)))
    logger.info(str(len(rslt_ef)))
    # ------ check slope
    equation = hull.equations
    new_equation = []
    for i in range(len(equation)):
        slope = -1*equation[i][0] / equation[i][1]
        if (slope != 0):
            new_equation.append([equation[i][0],equation[i][1],equation[i][2]])
    equation = new_equation
    logger.info('rate_x:'+str(rate_x))
    logger.info('rslt_ef:'+str(rslt_ef))

    # ------ remove end_point
    for i in range(2):
        del rate_x[0]
        del rslt_ef[0]

    # logger.info('after rate_x:'+str(rate_x))
    # logger.info('after rslt_ef:'+str(rslt_ef))
    # logger.info(str(len(rate_x)))
    # logger.info(str(len(rslt_ef)))

    # ------ calc hull distance(hdis)
    hdist = {}
    #for i in range(rin.tot_struc-c_struc_num,rin.tot_struc):
    for i in range(rin.tot_struc):
        # ------ check hull distance
        if(cal_hull_distance2d(rate_x[i],rslt_ef[i],equation)<0):
            raise SystemExit(1) #temporary
        
        hdist[i] = cal_hull_distance2d(rate_x[i],rslt_ef[i],equation)

    # ------ sort hdist
    logger.info('hdist:'+str(hdist))
    hdist = dict(sorted(hdist.items(),key=lambda x:x[1]))
    logger.info('len(hdist): '+str(len(hdist)))

    # ------ plot ID and save convex hull
    #plot ID(hdist == 0)
    for id in hdist.keys():
        if hdist[id] == 0.0:
            plt.text(rate_x[id]*0.95,rslt_ef[id]*1.55,f"ID:{id}",size='x-small')
            
    (gen, id_queueing,id_running) = pkl_data.load_ea_id()
    plt.savefig(f"./data/convex_hull_gen{gen}.png")

    # with open('./data/asc_hdist','a') as asc_h:
    #     asc_h.write(f"----- Gen{gen} -----\n")
    #     asc_h.write('ID'+'\t'+'hull distance'+'\n')
    #     asc_h.writelines("\n".join(str(i)+"\t"+str(j) for i,j in hdist.items()))
    #     asc_h.write('\n')

    return hdist