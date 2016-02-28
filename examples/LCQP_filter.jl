


function l1_filter2(y,trend=2,shift=0,spike=0,lambda=0,rho=0,gamma=0,l1l2=0,gap=0.03,time=1000)
    
    n=size(y)

    
    if (n[1]>1 && n[2]>1) 
         return :BoundError 
    else
        n = max(n[1],n[2])    
    end
    
    isempty(size(lambda)) == 1 ? lambda = lambda*ones(n,1) : true
    isempty(size(gamma))  == 1 ? gamma=gamma*ones(n,1)     : true
    isempty(size(rho))    == 1 ? rho=rho*ones(n,1)         : true


    #m = Model(solver = CbcSolver(ratioGap=gap, logLevel=0, seconds=time) );
    m = Model(solver = GurobiSolver() );
    
    #GurobiSolver(MIPGap = gap , TimeLimit=time, MIPFocus=1, SolutionLimit = 1)
    
    @defVar(m, l1_spike >= 0 )
    @defVar(m, l1_shift >= 0 )
    @defVar(m, l1_trend >= 0 )
    @defVar(m, u[t in 1:n] )
    @defVar(m, absU[t=1:n]>=0 )
    @defVar(m, w[t=1:n])
    @defVar(m, absW[t=2:n]>=0)

    @defVar(m, x[t=1:n])
    @defVar(m, absX[t=2:(n-1)]>=0)

    @defVar(m, 0 <= a <= 10)


    @defVar(m, goal[ t = 1:n])
    if l1l2 == 0
        @defVar(m, l1_goal[ t = 1:n]>=0)
    end
  
    for t in 1:n
        @addConstraint(m,  absU[t] >= u[t] )
        @addConstraint(m, -absU[t] <= u[t] )
    end
    @addConstraint(m, l1_spike == sum{ gamma[t]*absU[t] , t = 1:n})

    for t in 2:n
        @addConstraint(m,  absW[t] >= w[t]-w[t-1] )
        @addConstraint(m, -absW[t] <= w[t]-w[t-1] )
    end
    @addConstraint(m,  l1_shift == sum{ rho[t]*absW[t] , t=2:n } )

    for t in 2:(n-1)
        @addConstraint(m,  absX[t] >= x[t-1]-2*x[t]+x[t+1] )
        @addConstraint(m, -absX[t] <= x[t-1]-2*x[t]+x[t+1] )
    end
    @addConstraint(m, l1_trend == sum{ lambda[t]*absX[t] , t= 2:(n-1) } )

    for t in 1:n
        spike != 1 ? @addConstraint(m, 0 == u[t]) : true
        shift != 1 ? @addConstraint(m, 0 == w[t]) : true
        trend != 1 ? @addConstraint(m, 0 == x[t]) : true
    end
    trend != 2 ? @addConstraint(m, 0 == a ) : true

    for t in 2:(n-1)
        @defExpr(pregoal, y[t])
        spike == 1 ? @defExpr(pregoal, pregoal - u[t]) : true
        shift == 1 ? @defExpr(pregoal, pregoal - w[t]) : true
        trend == 1 ? @defExpr(pregoal, pregoal - x[t]) : true
        trend == 2 ? @defExpr(pregoal, pregoal - t*a ) : true
        @addConstraint(m, goal[t] == pregoal )
    end

    if l1l2==0
        for t = 1:n
            @addConstraint(m,  l1_goal[t] >= goal[t] )
            @addConstraint(m, -l1_goal[t] <= goal[t] )
        end
        @setObjective(m, Min, sum{l1_goal[t],t=2:(n-1)}      +l1_spike+l1_shift+l1_trend+lambda[1]*a )
    else
        @setObjective(m, Min, sum{goal[t]*goal[t],t=1:(n)} +l1_spike+l1_shift+l1_trend+lambda[1]*a )
    end

    status=solve(m)
    
    if status != :Optimal
        return status
    else 
        xx=zeros(n,1)
        ww=zeros(n,1)
        uu=zeros(n,1)
        aa=zeros(n,1)
        
        for t = 1:n
            xx[t]=getValue(x[t])
            ww[t]=getValue(w[t])
            uu[t]=getValue(u[t])
            aa[t]=t*getValue(a)
        end
        
       return xx, ww, uu, aa
    end
 
end#end func



function paramfree_l1(y,gridSize=100)
            xx_f=0
            uu_f=0
            ww_f=0
            aa_f=0
    n=size(y)
    
    (n[2]>n[1]) ? y= y' : true
    (n[1]>1 && n[2]>1) ? (return :BoundError) : n=max(n[1],n[2]) 
    
    grid=(logspace(0,1,gridSize)-1)/9
    
    Dtr=[eye(n-2) zeros(n-2,2)]-2*[zeros(n-2,1) eye(n-2) zeros(n-2,1)]+[zeros(n-2,2) eye(n-2) ]
    Dsh=[eye(n-1) zeros(n-1,1)]-1*[zeros(n-1,1) eye(n-1)]
    Dsp=eye(n)
    D=[Dtr,Dsh,Dsp]

    lambdaMax=maxabs((Dtr*Dtr')\Dtr*y)
    rhoMax=maxabs((Dsh*Dsh')\Dsh*y)*1
    gammaMax=maxabs((Dsp*Dsp')\Dsp*y)
    totalMax=maxabs((D*D')\D*y)
    
    bic_old = Inf
    
    #for l in lambdaMax*grid, r in rhoMax*grid
    for r in rhoMax*grid
        l=0
    #l1_filter2(y,trend=2,shift=0,spike=0,lambda=0,rho=0,gamma=0,l1l2=0,gap=0.03,time=1000)
        xx, ww, uu, aa = l1_filter2(y,0,1,0, l , r , gammaMax ,1 )
        epsilon=10^-5
        numParam = sum(abs(Dsp*(uu)).>epsilon) + sum(abs(Dsh*(ww)).>epsilon) + sum(abs(Dtr*(xx)).>epsilon)
        print(numParam)
       print("\n")
        SSE = sum((y-(xx + ww + uu + aa)).^2) 
         print(SSE)
        print("\n")
        bic = n*log(SSE)+numParam*log(n)
        
        if bic < bic_old
            
            xx_f=xx
            uu_f=uu
            ww_f=ww
            aa_f=aa
            bic_old = bic
            print(numParam)
            print("\n")
        end
        
    end
   
    return xx_f, ww_f, uu_f, aa_f

    
end





