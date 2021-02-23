#!/usr/bin/python3
from pylab import *
import pickle
from molecule import Molecule
import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
import matplotlib.patches as patches


PETRI_R = 10.0
R = 1.0

def vary_chi_plot():
    path = 'output/refuelling__fixed_env__motion__vary_chi'
    with open(path+'/data.p', 'rb') as handle:
        data = pickle.load(handle)
    ts = data['ts']
    ind = data['ind']

    ststs = data['trial_start_stops']#### THIS HSOULD BE trial_start_stops
    
    #chis = np.geomspace(0.5,20,21)
    chis = data['chis']
    # print(ind['x_h'])
    # quit()
    figure(figsize=(4.5,2.6))
    xs = []
    ys = []
    for trial_i,(start,stop) in enumerate(ststs):
        print(start,stop)
        x = ind['x_h'][stop-2]
        y = ind['y_h'][stop-2]
        ## show each trajectory
        # plot(ind['x_h'][start:stop],
        #      ind['y_h'][start:stop],'-',label=chis[trial_i],color=f'{float(trial_i)/len(ststs)}')
        # plot(x,y,'o')
        # title(f'{chis[trial_i]}')
        xs.append(x)
        ys.append(y)
        # plot(trial_i,
        #      y,'x-',label=chis[trial_i])
        #print(ind['x_h'][start:stop])
    # show()
    # final_its = np.array([stop-2 for _,stop in ststs],dtype=np.int)
    # print(final_its)
    # bar(chis,np.array(ind['x_h'])[final_its])
    plot(chis,xs,'ks-',ms=3,label='final x')
    #plot(chis,ys,'s-',label='final y')
    xscale('log')
    xlabel('$\chi_p$')
    ylabel('final $\\mathbf{q}_x$')
    title('final position of MOD')

    text(0.01,0.85,f'$[H] \\ll [P]$',transform=gca().transAxes,fontsize=7,verticalalignment='center',color='#0000ff')
    text(0.01,0.675,f'$[H] < [P]$',transform=gca().transAxes,fontsize=7,verticalalignment='center',color='#5500cc')
    text(0.01,0.5,f'$[H] = [P]$',transform=gca().transAxes,fontsize=7,verticalalignment='center',color='#bb00bb')
    text(0.01,0.325,f'$[H] > [P]$',transform=gca().transAxes,fontsize=7,verticalalignment='center',color='#cc0055')
    text(0.01,0.15,f'$[H] \\gg [P]$',transform=gca().transAxes,fontsize=7,verticalalignment='center',color='#ff0000')
    
    #legend()
    #show()
    tight_layout()
    savefig(f'plots/vary_chi.png',dpi=300)



def flow_plot():
    path = 'output/no_refuelling__fixed_env__motion'
    with open(path+'/data.p', 'rb') as handle:
        data = pickle.load(handle)
    ts = data['ts']
    ind = data['ind']

    figure('polar',figsize=(4,3))
    gs = gridspec.GridSpec(1,1)

    ## plot polar concentrations
    thetas = np.linspace(0,np.pi*2,48,endpoint=False)
    thetas = hstack([thetas,thetas[0]])
    ax = plt.subplot(gs[0,0], projection='polar')
    ax.xaxis.grid(False,color='none',linestyle='-')
    ax.yaxis.grid(False,color='none',linestyle='-')
    it = 50
    time = ts[it]
    for m in ['A','P','H'] :
        m_h = np.array(ind['c_hs'][m])
        scale = 1.0
        # if m == 'A':
        #     scale =
        #     # plot(ts[start_it:],scale*ys[start_it:],label=f'${scale}\\times[{m}]$',color=col)
        # else:
        #     scale = 0.0
        concs = m_h[it,:]
        concs = hstack([concs,concs[0]])
        concs = ((concs-np.min(concs))/(np.max(concs)-np.min(concs)))# + np.min(concs)#+0.01
        concs +=0.1
        concs *= 0.9
        #concs = log(scale*concs)
        plot(thetas,concs,'-',label=f'$[{m}]$',lw=2)
        ax.set_rticks([])  # less radial ticks
        ax.set_xticks(np.linspace(0,np.pi*2,4,endpoint=False))  # less radial ticks
        ax.set_xticklabels(['0','$\\pi/2$','$\\pi$','$3\\pi/2$'])  # less radial ticks
        #ax.set_rticklabels([])
        #ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line

        #plot(ts[it],ys[it],)
        # for r in range(np.shape(m_h)[1]) :
        #     plot(ts,m_h[:,r])

    for theta,flow in zip(thetas[::2],ind['marangoni_flux_h'][it][::2]) :
        #plot([t,t+0.1],[1.1,1.1])
        # r = 1.0
        # x = r*np.cos(t)
        # y = r*np.sin(t)
        # mag = 0.1
        # dx,dy = mag*np.cos(t+np.pi/2), mag*np.cos(t+np.pi/2),
        arrow(theta,1.1,flow*4000,0.0,head_width=0.04,fc='k')
        # a3 = patches.FancyArrowPatch((-0.4,-0.6), (0.4,-0.6))#,connectionstyle="arc3,rad=.5")
        # ax.add_patch(a3)
    ylim(0,1.2)

    text(-0.2,1.0,f't = {np.round(time)}',transform=ax.transAxes)
    legend(loc='right',bbox_to_anchor=(0.1, 0., 0.0, 0.0),fontsize=9)
    # cs = np.array([[section.concentrations[m]
    #                 for m in ind.model.surface_chemistry.molecules.values()]
    #                   for section in ind.sections])
    # ts = hstack([thetas,thetas[0]])
    # ts = vstack([ts,]*shape(cs)[1])
    # cs = vstack([cs[:,:],cs[0,:][newaxis]])
    # ax.plot(ts.T, log(cs))
    # ax.set_rmax(2)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
    #ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
    #ax.grid(True)
    tight_layout()
    #show()
    savefig(f'plots/flow.png',dpi=300)

def tricol_trajs():
    for exp_title in ['tricol_trajectories']:
        path = 'output/'+exp_title
        with open(path+'/data.p', 'rb') as handle:
            data = pickle.load(handle)
        ts = data['ts']
        ind = data['ind']
        ststs = data['traj_start_stops'] #### THIS HSOULD BE trial_start_stops
        # ststs = [ststs[i] for i in [7]]
        #ststs = [ststs[i] for i in [0,1,2,3,8]]
        #ststs = [ststs[i] for i in [4,5,6,9]]


        figure('spins')
        figure('vel_mag')
        mean_late_vels = []
        mean_late_dAs = []
        for trial_i,(start,stop) in enumerate(ststs):
            stop = stop-1
            dA = []
            vel_mags = []
            for seg_i in range(start,stop-1):
                x1 = ind['x_h'][seg_i]
                x2 = ind['x_h'][seg_i+1]
                y1 = ind['y_h'][seg_i]
                y2 = ind['y_h'][seg_i+1]
                a1 = np.arctan2(y1,x1)
                a2 = np.arctan2(y2,x2)
                diff = a2-a1
                if diff > np.pi  :
                    diff -= np.pi*2
                if diff < -np.pi :
                    diff += np.pi*2
                # diff = sign(diff)
                mean_x = (x1+x2)/2
                mean_y = (y1+y2)/2
                r = np.sqrt(mean_x**2 + mean_y**2)
                dA.append(diff)
                vel_mags.append(np.sqrt((y2-y1)**2 + (x2-x1)**2))
                #plot(diff)#[x1,x2],[y1,y2],alpha=0.1)
            figure('spins')
            plot(dA)
            figure('vel_mag')
            mean_late_vel = mean(vel_mags[-500:])
            mean_late_dA  = mean(dA[-500:])
            mean_late_vels.append(mean_late_vel)
            mean_late_dAs.append(mean_late_dA)
            plot(mean_late_dA,mean_late_vel,'o')
            yscale('log')
            xlabel('mean late dA')
            ylabel('mean late vel')
            #ylim(0,0.3)

        figure('spins')
        savefig(f'plots/{exp_title}_spins.png',dpi=300)
        figure('vel_mag')
        savefig(f'plots/{exp_title}_vel_mag.png',dpi=300)


        ## all traj plot
        figure(figsize=(4,4))
        ax = gca()
        for trial_i,(start,stop) in enumerate(ststs):
            stop = stop-1
            plot(ind['x_h'][start:stop-1],
                 ind['y_h'][start:stop-1],alpha=0.5)
            plot(ind['x_h'][stop-1],
                 ind['y_h'][stop-1],'ko',alpha=0.5)
            # plot(ind['x_h'][start],
            #      ind['y_h'][start],'kx',alpha=1.0)
            color = 'k'
            dA = mean_late_dAs[trial_i]
            if dA < 0.0 :
                color ='r'
            text(ind['x_h'][start],
                 ind['y_h'][start],f'{dA:0.2f}',alpha=1.0,fontsize=6,horizontalalignment='center',verticalalignment='center',color=color)

        petri = plt.Circle((0., 0.), PETRI_R,
                            color='None',ec='k', lw=4, clip_on=True)
        gca().add_artist(petri)
        k = 1.05*PETRI_R
        #k = 0.6*r
        ax.axis('off')
        tight_layout()
        xlim(-k,k)
        ylim(-k,k)
        savefig(f'plots/{exp_title}_scene.png',dpi=300)


        ## ##################### concentration plots
        for trial_i,(start,stop) in enumerate(ststs):
            stop = stop-1
            figure(figsize=(5.0,2.7))
            xlabel('t')

            for m in ['A','b','c','x','y','z'] :
                cfg = {
                    'A' : {'color':'c'},
                    'b' : {'color':'0.5'},
                    'c' : {'color':'m'},
                    'x' : {'color':'r'},
                    'y' : {'color':'g'},
                    'z' : {'color':'b'},
                }
                m_h = np.array(ind['c_hs'][m])
                ys = np.mean(m_h,axis=1)
                plot(ts[start:stop],ys[start:stop],label=f'[{m}]',**cfg[m])
                print('here')

                # ## plot all sections
                # for r in range(np.shape(m_h)[1]) :
                #     plot(ts,m_h[:,r])
            # ylim(0,0.00058)
            yscale('log')
            xlim(ts[start],ts[stop-1])
            legend(fontsize=9,loc='upper right',framealpha=1.0)
            savefig(f'plots/{exp_title}_timeseries_{trial_i}.png',dpi=300)

        quit()



def tricol_plots() :
    for exp_title in ['tricol']:
        path = 'output/'+exp_title
        with open(path+'/data.p', 'rb') as handle:
            data = pickle.load(handle)
        ts = data['ts']
        ind = data['ind']

        #############################################33
        ### Scene plots
        env = data['env']
        fig = figure(figsize=(3.2,3.2))
        ax = gca()
        r = PETRI_R
        c = env['x']
        #imshow(c,cmap='gray',extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')
        #imshow(c,cmap='gray',extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')
        res = np.shape(c)[0]
        rgb = np.zeros((res,res,3))
        rgb[:,:,0] = env['x']/np.max(env['x'])
        rgb[:,:,1] = env['y']/np.max(env['y'])
        rgb[:,:,2] = env['z']/np.max(env['z'])
        # k = 20000.0
        # rgb[:,:,0] = env['H']*k
        # rgb[:,:,1] = env['P']*k
        # rgb[:,:,2] = 2*env['P']*k
        rgb[rgb[:,:,:] > 0.9] = 0.9
        rgb[env['mask'],:] = 1.0

        imshow(rgb,extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')

        petri = plt.Circle((0., 0.), r,
                            color='None',ec='k', lw=4, clip_on=True)
        ax.add_artist(petri)

        ax.axis('off')
        tight_layout()
        k = 1.05*r
        #k = 0.6*r
        xlim(-k,k)
        ylim(-k,k)
        plot(ind['x_h'][1:],ind['y_h'][1:],'w--',lw=1.)

        mod = plt.Circle((ind['x_h'][-1], ind['y_h'][-1]), R,
                         color='0.6', ec='k', clip_on=False, zorder=10)
        ax.add_artist(mod)

        #title(f'time = {ts[-1]:0.2}')
        title(f'time = {np.round(ts[-1])}')
        savefig(f'plots/{exp_title}_scene.png',dpi=300)


        ## ##################### concentration plots
        figure(figsize=(5.0,2.7))
        # plot(ts,[peak_y-y for y in ind['y_h']],'k',lw=2.5,label='y')
        # ylim(-2,12.0)
        # ylabel('y')
        xlabel('t')


        start_it = 2
        # ax = twinx()
        for m in ['A','b','c','x','y','z'] :
            cfg = {
                'A' : {'color':'c'},
                'b' : {'color':'0.5'},
                'c' : {'color':'m'},
                'x' : {'color':'r'},
                'y' : {'color':'g'},
                'z' : {'color':'b'},
            }
            m_h = np.array(ind['c_hs'][m])
            ys = np.mean(m_h,axis=1)
            plot(ts[start_it:],ys[start_it:],label=f'[{m}]',**cfg[m])
            # ## plot all sections
            # for r in range(np.shape(m_h)[1]) :
            #     plot(ts,m_h[:,r])
        # ylim(0,0.00058)
        xlim(0,75)
        legend(fontsize=9,loc='upper right',framealpha=1.0)
        savefig(f'plots/{exp_title}_timeseries.png',dpi=300)


def recreate_plots() :
    for exp_i,(exp_title,label,letter) in \
        enumerate([['no_refuelling__fixed_env__motion','No Refuelling (Static Env.)','A'],
                   ['no_refuelling__fixed_env__no_motion','No Refuelling (No Motion, Static Env.)','B'],
                   ['no_refuelling__dynamic_env__no_motion','No Refuelling (No Motion, Dyn. Env.)','C'],
                   ['no_refuelling__dynamic_env__motion','No Refuelling (Dyn. Env.)','D'],
                   ####
                   ['refuelling__fixed_env__motion','Refuelling (Static Env.)','A'],
                   ['refuelling__fixed_env__no_motion','Refuelling (No Motion, Static Env.)','B'],
                   ['refuelling__dynamic_env__no_motion','Refuelling (No Motion, Dyn. Env.)','C'],
                   ['refuelling__dynamic_env__motion','Refuelling (Dyn. Env.)','D']]):
        path = 'output/'+exp_title
        with open(path+'/data.p', 'rb') as handle:
            data = pickle.load(handle)

        ts = data['ts']
        ind = data['ind']
        peak_y = 0.0

        final_x = ind['x_h'][-1]
        final_y = ind['y_h'][-1]
        print(f'final resting point for "{exp_title}" is ({final_x},{final_y})')

        start_it = 2

        figure(figsize=(5.0,2.7))
        # plot(ts,[peak_y-y for y in ind['y_h']],'k',lw=2.5,label='y')
        # ylim(-2,12.0)
        # ylabel('y')
        xlabel('t')

        ylabel('$[P],[H],[A],\Delta \Gamma_{max}$')
        for m in ['A','P','H'] :
            m_h = np.array(ind['c_hs'][m])
            ys = np.mean(m_h,axis=1)
            # if m == 'A' :
            #     plot([0,0],[1E-5,1E-5],label=f'[{m}]')
            # else :
            plot(ts[start_it:],ys[start_it:],label=f'[{m}]')
            # ## plot all sections
            # for r in range(np.shape(m_h)[1]) :
            #     plot(ts,m_h[:,r])

        # # plot 'H' differently because of scaling
        # scale = 1.0#0.001
        # h_h = np.array(ind['c_hs']['H'])
        # plot(ts,np.mean(h_h,axis=1)*scale,label=f'$[h]\\times{scale}$')

        plot(ts[start_it:],[y for y in ind['st_diff_h']][start_it:],
             'k--',lw=1,label='$\Delta \Gamma_{max}$')
        plot([0,0],[1E-15,1E-15],'s',color='0.5',label=f'$d[A]/dt$')

        ax = gca()
        text(0.985,1.02,label,horizontalalignment='right',
             verticalalignment='bottom',transform=ax.transAxes)

        text(-0.06,1.12,letter,horizontalalignment='center',
             verticalalignment='top',transform=ax.transAxes,fontsize=20, weight='bold')

        xlim(0,ts[-1])
        ylim(1E-8,2E-1)
        yscale('log')
        xlabel('t')
        if letter == 'A' :
            legend(fontsize=9,loc='lower right')

        # ax = twinx()
        # ylabel('$d[A]/dt$')
        for m in ['A'] :
            m_h = np.array(ind['c_hs'][m])
            ys = np.mean(m_h,axis=1)
            # plot(ts[start_it:],ys[start_it:],label=f'[{m}]')
            #plot(ts[start_it:],np.gradient(ys[start_it:])/0.05,'b:') ## ciao
            y = np.array(np.gradient(ys[start_it:])/0.05)
            print(f'last dadt for {exp_title}: {y[-1]}')
            fill_between(ts[start_it:],0*y,y2=y,color='0.5',edgecolor='0.5',alpha=0.2) ## ciao
            
            # ## plot all sections
            # for r in range(np.shape(m_h)[1]) :
            #     plot(ts,m_h[:,r])
        #ylim(0.,0.0004)
        yscale('log')
        tight_layout()
        savefig(f'plots/{exp_title}_timeseries.png',dpi=300)


        #############################################
        ### Scene plots
        env = data['env']
        fig = figure(figsize=(3.2,3.2))
        ax = gca()
        r = PETRI_R
        c = env['H']
        #imshow(c,cmap='gray',extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')
        #imshow(c,cmap='gray',extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')
        res = np.shape(c)[0]
        rgb = np.zeros((res,res,3))
        rgb[:,:,0] = env['H']/np.max(env['H'])
        rgb[:,:,1] = env['P']/np.max(env['P'])
        rgb[:,:,2] = 2*env['P']/np.max(env['P'])
        # k = 20000.0
        # rgb[:,:,0] = env['H']*k
        # rgb[:,:,1] = env['P']*k
        # rgb[:,:,2] = 2*env['P']*k
        #rgb[rgb[:,:,:] > 0.9] = 0.9
        rgb[env['mask'],:] = 1.0

        imshow(rgb,extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')

        petri = plt.Circle((0., 0.), r,
                            color='None',ec='k', lw=10, clip_on=True)
        ax.add_artist(petri)

        #ax.axis('off')
        tight_layout()
        #k = 1.05*r
        if exp_i < 4:
            y = -5.1
            x = -2.5
            r = 3.6
            xlim(x,x+r)
            ylim(y,y+r)
        else :
            k = 0.6*r
            xlim(-k,k)
            ylim(-k,k)
        plot(ind['x_h'][10:],ind['y_h'][10:],'w--',lw=1.5)

        mod = plt.Circle((ind['x_h'][-1], ind['y_h'][-1]), R,
                         color='0.6', ec='k')#, clip_on=False, zorder=10)
        ax.add_artist(mod)

        #title(f'time = {ts[-1]:0.2}')
        title(f'time = {np.round(ts[-1])}')
        savefig(f'plots/{exp_title}_scene.png',dpi=300)

        #### chemical distribution plots
        figure(figsize=(8,4))
        #get_current_fig_manager().window.wm_geometry("-0+0")
        n = len(ind['c_hs'].keys())
        gs = gridspec.GridSpec(n,1)

        for c_h_i,key in enumerate(ind['c_hs'].keys()) :
            ax = plt.subplot(gs[c_h_i,0])
            ch = np.array(ind['c_hs'][key])
            imshow(ch.T,aspect='auto',extent=(0,ts[-1],0,np.pi*2),origin='lower',cmap='gray_r')
            # ax.set_title(str(key))
            ylabel(str(key),rotation=0,fontsize=14)
            if c_h_i != 2:
                xticks([])
            # else :
            #     ylabel('t')
            yticks([0,np.pi,2.0*np.pi],fontsize=12)
            gca().set_yticklabels(['0','$\\pi$','$2\pi$'])
            #xlim(0,40)
            if c_h_i == len(ind['c_hs'].keys())-1:
                xlabel('t',fontsize=14)
        tight_layout()

        savefig(f'plots/{exp_title}_surface.png',dpi=300)


    # show()


def recreate_where_a_made() :
    for exp_title,label,letter in \
    zip(['refuelling__fixed_env__motion'],#'refuelling'],
        ['',''],
        ['','']):
        path = 'output/'+exp_title
        with open(path+'/data.p', 'rb') as handle:
            data = pickle.load(handle)

        ts = data['ts']
        ind = data['ind']
        peak_y = 0.0

        start_it = 2

        # ################# All SECTIONS
        # figure(figsize=(8.0,5.0))
        # xlabel('t')
        # # ax = twinx()
        # for m,col in zip(['A','P','H'],
        #                  ['r','g','b']):
        #     m_h = np.array(ind['c_hs'][m])
        #     ys = np.mean(m_h,axis=1)
        #     scale = 1.0
        #     if m == 'A':
        #         scale = 0.05
        #     #plot(ts[start_it:],scale*ys[start_it:],label=f'${scale}\\times[{m}]$',color=col)
        #     ## plot all sections
        #     for r in range(np.shape(m_h)[1]) :
        #         plot(ts,m_h[:,r]*scale,color=col,alpha=0.2)

        ################# ONLY SELECTED SECTIONS
        figure(figsize=(8.0,5.0))
        xlabel('t')
        # ax = twinx()
        # for m,col in zip(['P','H'],
        #                  ['g','b']):
        #     m_h = np.array(ind['c_hs'][m])
        #     ys = np.mean(m_h,axis=1)
        #     scale = 1.0
        #     if m == 'A':
        #         scale = 0.05
        #     #plot(ts[start_it:],scale*ys[start_it:],label=f'${scale}\\times[{m}]$',color=col)
        #     ## plot all sections
        #     for r in [12,36] :
        #         plot(ts,m_h[:,r]*scale,color=col,alpha=0.2)
        # xlim(8,15)
        # # # plot 'H' differently because of scaling
        # # scale = 1.0#0.001
        # # h_h = np.array(ind['c_hs']['H'])
        # # plot(ts,np.mean(h_h,axis=1)*scale,label=f'$[h]\\times{scale}$')

        # plot(ts[start_it:],[y for y in ind['st_diff_h']][start_it:],
        #      'k--',lw=1,label='$\Delta \Gamma_{max}$')
        yscale('log')
        legend()

        savefig(f'plots/{exp_title}_detail.png',dpi=300)

        ### where is a being made?
        figure(figsize=(5.5,3.25))
        #get_current_fig_manager().window.wm_geometry("-0+0")
        p = np.array(ind['c_hs']['P']).T
        h = np.array(ind['c_hs']['H']).T
        # a = np.array(ind['c_hs']['A']).T
        # for col in range(np.shape(a)[1]) :
        #     a[:,col] /= np.max(a[:,col])
        rgb = np.zeros((shape(p)[0],shape(p)[1],3))
        # rgb[:,:,0] = p / np.max(p)
        # rgb[:,:,1] = h / np.max(h)
        rgb[:,:,0] = h*p / np.max(h*p)
        rgb[:,:,1] = h*p / np.max(h*p)
        rgb[:,:,2] = h*p / np.max(h*p)

        #imshow(rgb,aspect='auto',extent=(0,ts[-1],0,np.pi*2),origin='lower',cmap='gray_r')
        imshow(h*p / np.max(h*p),
               aspect='auto',extent=(0,ts[-1],0,np.pi*2),origin='lower',cmap='gray_r')
        # imshow(a,
        #        aspect='auto',extent=(0,ts[-1],0,np.pi*2),origin='lower',cmap='gray')
        #colorbar()
        plot(ts,np.argmax(rgb[:,:,2],axis=0)/48.0*np.pi*2.0,'.',color='#ffff00',ms=5,
             label='$max_{t=t} [P][H]$')

        # xmax = 33.0
        # xmin = 42.0
        # dt = data['DT']
        # xmax_it = int(xmax/dt)
        # xmin_it = int(xmin/dt)

        ys = ind['y_h']
        xs = ind['x_h']
        dys = np.gradient(ys)
        dxs = np.gradient(xs)
        vel_data = np.array([(np.arctan2(y,x)+2.0*np.pi)%(np.pi*2.0) for y,x in zip(dys,dxs)])
        # abs_d_data = np.abs(np.diff(vel_data))
        # mask = np.hstack([ abs_d_data > abs_d_data.mean()+3*abs_d_data.std(), [False]])
        # masked_data = np.ma.MaskedArray(vel_data, mask)
        # plot(ts,masked_data,'c-')
        plot(ts,vel_data,'s',color='#00ffff',ms=2.5,label='tan$^{-1}\\left (\\frac{d\\mathbf{q}}{dt}\\right )$')

        # figure()
        # plot(ts,xs)
        # show()
        # quit()
        
        xmax_it = np.argmax(xs)
        xmin_it = np.argmin(xs)
        print(f'MAX: at it {xmax_it}, xs={xs[xmax_it]}')
        print(f'MIN: at it {xmin_it}, xs={xs[xmin_it]}')
        yticks([0,
                np.pi/2,
                np.pi,
                3.0*np.pi/2.0],fontsize=12)
        gca().set_yticklabels(['$0 (\\rightarrow)$','$\\pi/2 (\\uparrow)$','$\\pi (\\leftarrow)$','$3\\pi/2 (\\downarrow)$'])#,'$2\pi$'])
        #gca().set_yticklabels(['$0$','$\\pi/2$','$\\pi$','$3\\pi/2$'])#,'$2\pi$'])
        ylabel('$\\theta$')
        xlabel('t')
        #xlim(0,25)
        l = legend(bbox_to_anchor=(0.0, 0., 1.0, 1.5))
        frame = l.get_frame()
        frame.set_facecolor('0.66')
        frame.set_edgecolor('0.2')

        title(r'$[P][H] \propto \frac{d[A]}{dt}$')
        tight_layout()

        savefig(f'plots/{exp_title}_where_a_made.png',dpi=300)

        figure()
        it_with_x_maximum = xmax_it#130
        it_with_x_minimum = xmin_it#190
        t_with_x_maximum = ts[it_with_x_maximum]
        t_with_x_minimum = ts[it_with_x_minimum]

        r_loc = 0
        l_loc = 24

        gs = gridspec.GridSpec(4,1)
        def plot_product(theta) :
            plot(ts,p[theta,:],label='[P]')
            plot(ts,h[theta,:],label='[H]')
            ax = twinx()
            plot(ts,[x*y for x,y in zip(p[theta,:],h[theta,:])],label='[P][H]',color='k')
            #xlim(7,20)
        ax = plt.subplot(gs[0,0])
        plot_product(l_loc)
        ax = plt.subplot(gs[1,0])
        plot_product(r_loc)
        ax = plt.subplot(gs[2,0])
        plot(ts,ind['x_h'])
        plot(ts,ind['y_h'])
        ax = plt.subplot(gs[3,0])
        title('Use this one to show Correlation between p-production and direction travelled. Green is higher than black when moving up.')
        plot(ts,[x*y for x,y in zip(p[l_loc,:],h[l_loc,:])],label='[P][H]',color='g')
        plot(ts,[x*y for x,y in zip(p[r_loc,:],h[r_loc,:])],label='[P][H]',color='k')
        plot(ts,ind['x_h'])
        tight_layout()
        savefig(f'plots/{exp_title}_where_a_made_exp.png',dpi=300)

        def barchart(title,iteration):
            figure(figsize=(5.5,2.5))
            gs = gridspec.GridSpec(1,2,width_ratios=[0.6,0.4])
            ax = plt.subplot(gs[0,0])
            gap = 0.5
            xs = [1,2,3+gap,4+gap]
            ys = []
            labels = []

            ys.append(p[l_loc][iteration])
            labels.append('$[P]_{\\theta=\\pi}$')
            ys.append(p[r_loc][iteration])
            labels.append('$[P]_{\\theta=0}$')

            ys.append(h[l_loc][iteration])
            labels.append('$[H]_{\\theta=\\pi}$')
            ys.append(h[r_loc][iteration])
            labels.append('$[H]_{\\theta=0}$')

            xticks(xs)
            gca().set_xticklabels(labels)
            suptitle(f'At $q_x$-{title} ($t\\approx{np.round(ts[iteration],1)}$)                ',fontsize=10)
            bar(xs,ys,color=['#ffa500','#ffa500','g','g'])
            ylim(0.0,0.0001)

            ax = plt.subplot(gs[0,1])
            xs = [1,2]
            ys = []
            xlim(0.4,2.6)
            labels = []
            ys.append(p[l_loc][iteration]*h[l_loc][iteration])
            #labels.append('$\\frac{d[A]}{dt}_{\\theta=\\pi}$')
            labels.append('$[P][H]_{\\theta=\\pi}$')
            ys.append(p[r_loc][iteration]*h[r_loc][iteration])
            #labels.append('$\\frac{d[A]}{dt}_{\\theta=0}$')
            labels.append('$[P][H]_{\\theta=0}$')

            xticks(xs)
            gca().set_xticklabels(labels)
            bar(xs,ys,color='k')
            ylim(0,3E-9)
            tight_layout()

            savefig(f'plots/bar_{title}.png')

        barchart('maximum',it_with_x_maximum)
        barchart('minimum',it_with_x_minimum)
        
        print()
        print(f'At t_max_x ({t_with_x_maximum}), MOD is close to P-peak')
        print(f'left location')
        print(f'p[{l_loc}] = {p[l_loc][it_with_x_maximum]}')
        print(f'h[{l_loc}] = {h[l_loc][it_with_x_maximum]}')
        print(f'product = {p[l_loc][it_with_x_maximum]*h[l_loc][it_with_x_maximum]}')
        print(f'right location')
        print(f'p[{r_loc}] = {p[r_loc][it_with_x_maximum]}')
        print(f'h[{r_loc}] = {h[r_loc][it_with_x_maximum]}')
        print(f'product = {p[r_loc][it_with_x_maximum]*h[r_loc][it_with_x_maximum]}')


        print()
        print(f'At t_min_x ({t_with_x_minimum}), MOD is close to P-peak')
        print(f'left location')
        print(f'p[{l_loc}] = {p[l_loc][it_with_x_minimum]}')
        print(f'h[{l_loc}] = {h[l_loc][it_with_x_minimum]}')
        print(f'product = {p[l_loc][it_with_x_minimum]*h[l_loc][it_with_x_minimum]}')
        print(f'right location')
        print(f'p[{r_loc}] = {p[r_loc][it_with_x_minimum]}')
        print(f'h[{r_loc}] = {h[r_loc][it_with_x_minimum]}')
        print(f'product = {p[r_loc][it_with_x_minimum]*h[r_loc][it_with_x_minimum]}')
        #show()



def total_a_produced() :
    experiments = \
        [['no_refuelling__fixed_env__motion','No Refuelling (Static Env.)','1A'],
         ['no_refuelling__fixed_env__no_motion','No Refuelling (No Motion, Static Env.)','1B'],
         ['no_refuelling__dynamic_env__no_motion','No Refuelling (No Motion, Dyn. Env.)','1C'],
         ['no_refuelling__dynamic_env__motion','No Refuelling (Dyn. Env.)','1D'],
         ####
         ['refuelling__fixed_env__motion','Refuelling (Static Env.)','2A'],
         ['refuelling__fixed_env__no_motion','Refuelling (No Motion, Static Env.)','2B'],
         ['refuelling__dynamic_env__no_motion','Refuelling (No Motion, Dyn. Env.)','2C'],
         ['refuelling__dynamic_env__motion','Refuelling (Dyn. Env.)','2D']]

    figure(figsize=(9.0,3.0))
    
    results = dict()
    final_a = dict()
    for exp_i,(exp_title,label,letter) in enumerate(experiments) :
        path = 'output/'+exp_title
        with open(path+'/data.p', 'rb') as handle:
            data = pickle.load(handle)

        ts = data['ts']
        ind = data['ind']
        peak_y = 0.0

        start_it = 2
        end_it = 800
        
        h_h = np.sum( np.array(ind['c_hs']['H'])[start_it:end_it],axis=1)
        p_h = np.sum( np.array(ind['c_hs']['P'])[start_it:end_it],axis=1)
        a_h = np.sum( np.array(ind['c_hs']['A'])[start_it:end_it],axis=1)
        a_end = np.sum( np.array(ind['c_hs']['A'])[end_it-1])

        
        print(len(h_h))
        prod = h_h*p_h
        results[exp_title] = prod
        final_a[exp_title] = a_end
        linestyle = '-'
        if letter[0] == '2' :
            linestyle = '--'
        plot(ts[start_it:end_it],prod,label=f'{letter}: {label}',ls=linestyle)
        
    legend(loc='upper left',bbox_to_anchor=(1.0, 1.022,0.0,0.0),fontsize=9)#,bbox_transform=plt.gcf().transFigure)
    ylabel(r'$[P][H] \propto \frac{d[A]}{dt}$')
    xlabel('t')
    xlim(0.0,40.0)
    title(r'$A$ production')
    tight_layout()
    savefig(f'plots/a_production.png',dpi=300)

    figure(figsize=(5.0,2.7))
    labels = []
    final_as = []
    prods = []
    xs = []
    for exp_i,(exp_title,label,letter) in enumerate(experiments) :
        xs.append(exp_i)
        labels.append(letter)
        final_as.append(final_a[exp_title])
        prods.append(np.sum(results[exp_title]))
    print(xs)
    #bar(xs,prods)
    bar(xs,final_as)
    xticks(xs)
    xlabel('experiment')
    ylabel('$[A]_{t=40}$')
    title('Total A produced at $t=40$')
    # text(0.05,0.9,f'$t \\approx {ts[end_it]:0.1f}$',
    #      transform=gca().transAxes,fontsize=7,verticalalignment='center',color='k')
    gca().set_xticklabels(labels)
    tight_layout()
    savefig(f'plots/total_a_produced.png',dpi=300)

    # def barchart(title,iteration):
    #     figure(figsize=(5.5,2.5))
    #     gs = gridspec.GridSpec(1,2,width_ratios=[0.6,0.4])
    #     ax = plt.subplot(gs[0,0])
    #     gap = 0.5
    #     xs = [1,2,3+gap,4+gap]
    #     ys = []
    #     labels = []

    #     ys.append(p[l_loc][iteration])
    #     labels.append('$[P]_{\\theta=\\pi}$')
    #     ys.append(p[r_loc][iteration])
    #     labels.append('$[P]_{\\theta=0}$')

    #     ys.append(h[l_loc][iteration])
    #     labels.append('$[H]_{\\theta=\\pi}$')
    #     ys.append(h[r_loc][iteration])
    #     labels.append('$[H]_{\\theta=0}$')

    #     xticks(xs)
    #     gca().set_xticklabels(labels)
    #     suptitle(f'At $q_x$-{title} ($t\\approx{np.round(ts[iteration],1)}$)                ',fontsize=10)
    #     bar(xs,ys,color=['#ffa500','#ffa500','g','g'])
    #     ylim(0.0,0.0001)

    #     ax = plt.subplot(gs[0,1])
    #     xs = [1,2]
    #     ys = []
    #     xlim(0.4,2.6)
    #     labels = []
    #     ys.append(p[l_loc][iteration]*h[l_loc][iteration])
    #     #labels.append('$\\frac{d[A]}{dt}_{\\theta=\\pi}$')
    #     labels.append('$[P][H]_{\\theta=\\pi}$')
    #     ys.append(p[r_loc][iteration]*h[r_loc][iteration])
    #     #labels.append('$\\frac{d[A]}{dt}_{\\theta=0}$')
    #     labels.append('$[P][H]_{\\theta=0}$')

    #     xticks(xs)
    #     gca().set_xticklabels(labels)
    #     bar(xs,ys,color='k')
    #     ylim(0,3E-9)
    #     tight_layout()

    #     savefig(f'plots/bar_{title}.png')

    # barchart('maximum',it_with_x_maximum)
    # barchart('minimum',it_with_x_minimum)


#recreate_plots()
#recreate_where_a_made()
#flow_plot()
vary_chi_plot()

#total_a_produced()


