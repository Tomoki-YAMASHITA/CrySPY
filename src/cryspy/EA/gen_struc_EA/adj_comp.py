from logging import getLogger
import random
import numpy as np
from pymatgen.transformations.site_transformations import ReplaceSiteSpeciesTransformation
from ...IO import pkl_data
from ...IO import read_input as rin
from ...util.struc_util import get_nat


logger = getLogger('cryspy')


def convex_hull_check():
    '''
    convex hullのどの区間にどれくらいの構造があるかをカウントする
    とりあえずconvex hullを5分割して考える
    ''' 
    # どの区間にどれだけの構造があるかをカウントする
    nat_data, ratio_data = pkl_data.load_ea_vc_data()
    section = [0,0,0,0,0]
    for i in range(len(nat_data)):
        tmp_arr = ratio_data[i][0]
        #logger.info(f'ratio: {tmp_arr}')
        if 0 <= tmp_arr <0.2:
            section[0] +=1
        elif 0.2 <= tmp_arr <0.4:
            section[1] +=1
        elif 0.4 <= tmp_arr <0.6:
            section[2] +=1
        elif 0.6 <= tmp_arr <0.8:
            section[3] +=1
        elif 0.8 <= tmp_arr <=1.0:
            section[4] +=1
        logger.info(f'section: {section}')

    
    if rin.target == 'depop':
        tmp_arr = []
        for i in range(len(section)):
            if i == 0:
                tmp_arr.append(0)
            else:
                if section[tmp_arr[0]] > section[i]:
                    tmp_arr = []
                    tmp_arr.append(i)
                elif section[tmp_arr[0]] == section[i]:
                    tmp_arr.append(i)
        
        if len(tmp_arr) > 1:
            rslt = random.choice(tmp_arr)
        else:
            rslt = tmp_arr[0]
        logger.info(f'tmp_arr: {tmp_arr}')
        logger.info(f'select_section: {rslt}')

    # ---------- 構造が多いところに子構造を生成する場合の処理    
    elif rin.target == 'overpop':
        tmp_arr = []
        for i in range(len(section)):
            if i == 0:
                tmp_arr.append(i)
            else:
                # logger.info(f'i: {i}')
                # logger.info(f'{section[tmp_arr[0]]} < {section[i]}')
                if section[tmp_arr[0]] < section[i]:
                    tmp_arr = []
                    tmp_arr.append(i)
                elif section[tmp_arr[0]] == section[i]:
                    tmp_arr.append(i)
        if len(tmp_arr) > 1:
            rslt = random.choice(tmp_arr)
        else:
            rslt = tmp_arr[0]
        logger.info(f'tmp_arr: {tmp_arr}')
        logger.info(f'select_section: {rslt}')

    return rslt

def operation_atoms(method,child, section):
    """
    構造が少ないところをターゲットにして子構造生成を行う
    methondにはaddition, elimination, substitutionのいずれかを指定する

    section = 0 : 0<=x<0.2
    section = 1 : 0.2<=x<0.4
    section = 2 : 0.4<=x<0.6
    section = 3 : 0.6<=x<0.8
    section = 4 : 0.8<=x<=1.0
    """
    # どの区間に対して子構造を生成するかを決める
    # 候補が多い場合はランダムで一つ決める
    #section = convex_hull_check()
    
    nat,ratio = get_nat(child,rin.atype)
    ratio = ratio[0]
    # ---------- addition
    if method == 'addition':
        logger.info(' ---------- addition')
        cnt = 0
        logger.info(f'section: {section}')
        while True:
            coords = np.random.rand(3)
            if section == 0:
                if 0<= ratio < 0.2: 
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.2:
                    child.append(species=rin.atype[1], coords=coords)
                    
            if section == 1:
                if 0.2<= ratio < 0.4:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.4:
                    child.append(species=rin.atype[1], coords=coords)
                    
                else:
                    child.append(species=rin.atype[0], coords=coords)

            if section == 2:
                if 0.4<= ratio < 0.6:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.6:
                    child.append(species=rin.atype[1], coords=coords)
                    
                else:
                    child.append(species=rin.atype[0], coords=coords)
                
            if section == 3:
                if 0.6<= ratio < 0.8:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.8:
                    child.append(species=rin.atype[1], coords=coords)
                    
                else:
                    child.append(species=rin.atype[0], coords=coords)
            
            if section == 4:
                if 0.8<= ratio <= 1.0:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 1.0:
                    child.append(species=rin.atype[1], coords=coords)
                    
                else:
                    child.append(species=rin.atype[0], coords=coords)

            nat,ratio = get_nat(child,rin.atype)
            ratio = ratio[0]
            cnt += 1
            if cnt == rin.maxcnt_ea:
                #logger.info(f'return None')
                return None


    # ---------- elimination
    elif method == 'elimination':
        logger.info(' ---------- elimination')
        cnt = 0
        logger.info(f'section: {section}')
        # ------ prepare index for each atom type
        indx_each_type = []
        for a in rin.atype:
            indx_each_type.append(
                [i for i, site in enumerate(child)
                    if site.species_string == a])
        while True:
            
            if section == 0:
                if 0<= ratio < 0.2: 
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.2:
                    indx = indx_each_type[0]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                    
            if section == 1:
                if 0.2<= ratio < 0.4:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.4:
                    indx = indx_each_type[0]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                    
                else:
                    indx = indx_each_type[1]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)

            if section == 2:
                if 0.4<= ratio < 0.6:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.6:
                    indx = indx_each_type[0]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                    
                else:
                    indx = indx_each_type[1]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                
            if section == 3:
                if 0.6<= ratio < 0.8:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.8:
                    indx = indx_each_type[0]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                    
                else:
                    indx = indx_each_type[1]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
            
            if section == 4:
                if 0.8<= ratio <= 1.0:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 1.0:
                    indx = indx_each_type[0]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
                    
                else:
                    indx = indx_each_type[1]
                    rm = np.random.choice(indx,1,replace=False)
                    child.remove_sites(rm)
            from ..ea_child import check_vcnat
            nat_success = check_vcnat(child)
            if nat_success == False:
                return None
            nat,ratio = get_nat(child,rin.atype)
            ratio = ratio[0]
            cnt += 1
            if cnt == rin.maxcnt_ea:
                logger.info(f'return None')
                return None
        
    # ---------- substitution
    else:
        logger.info(' ---------- subsitution')
        cnt = 0
        logger.info(f'section: {section}')
        indx_each_type = []
        for a in rin.atype:
            indx_each_type.append(
                [i for i, site in enumerate(child)
                    if site.species_string == a])
            

        while True:
            if section == 0:
                if 0<= ratio < 0.2: 
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.2:
                    indx = indx_each_type[0]
                    rep = np.random.choice(indx,1,replace=False)
                    logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[1]})
                    child = transform.apply_transformation(child)
                    
            if section == 1:
                if 0.2<= ratio < 0.4:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.4:
                    indx = indx_each_type[0]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[1]})
                    child = transform.apply_transformation(child)
                    
                else:
                    indx = indx_each_type[1]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[0]})
                    child = transform.apply_transformation(child)

            if section == 2:
                if 0.4<= ratio < 0.6:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.6:
                    indx = indx_each_type[0]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[1]})
                    child = transform.apply_transformation(child)
                    
                else:
                    indx = indx_each_type[1]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[0]})
                    child = transform.apply_transformation(child)
                
            if section == 3:
                if 0.6<= ratio < 0.8:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 0.8:
                    indx = indx_each_type[0]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[1]})
                    child = transform.apply_transformation(child)
                    
                else:
                    indx = indx_each_type[1]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[0]})
                    child = transform.apply_transformation(child)
            
            if section == 4:
                if 0.8<= ratio <= 1.0:
                    if cnt == 0:
                        return None
                    else:
                        return child
                elif ratio >= 1.0:
                    indx = indx_each_type[0]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[1]})
                    child = transform.apply_transformation(child)
                    
                else:
                    indx = indx_each_type[1]
                    rep = np.random.choice(indx,1,replace=False)
                    #logger.info(f'rep: {rep}')
                    transform = ReplaceSiteSpeciesTransformation({rep[0]: rin.atype[0]})
                    child = transform.apply_transformation(child)

            from ..ea_child import check_vcnat
            nat_success = check_vcnat(child)
            if nat_success == False:
                return None
            nat,ratio = get_nat(child,rin.atype)
            ratio = ratio[0]
            cnt += 1
            if cnt == rin.maxcnt_ea:
                #logger.info(f'return None')
                return None








        
    