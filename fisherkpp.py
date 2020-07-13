import sys
from ngsolve import *
from ngsolve.solvers import *
from prodmesh import TensorProdMesh,CartSquare
from trefftzngs import *

def SolFisherKPP(mesh,order,uin,eps,steps,nbndc=[],nbnd="",acoeff=1,Dcoeff=1):
    N=uin.dim
    D=mesh.dim
    entropy = []
    V = H1(mesh,order=order)
    X = FESpace(([V for i in range(N)]))
    trial = X.TrialFunction()
    test = X.TestFunction()
    w = CoefficientFunction( tuple([trial[i] for i in range(N)]) )
    v = CoefficientFunction( tuple([ test[i] for i in range(N)]) )

    dtw = CoefficientFunction(tuple([grad(trial[i])[D-1] for i in range(N)]))
    dtv = CoefficientFunction(tuple([grad( test[i])[D-1] for i in range(N)]))
    dxw = [0]*(D-1)
    dxv = [0]*(D-1)
    for i in range((D-1)):
        dxw[i] = CoefficientFunction(tuple([grad(trial[j])[i] for j in range(N)]))
        dxv[i] = CoefficientFunction(tuple([grad( test[j])[i] for j in range(N)]))

    u = (1/(1+sum(exp(w))))*exp(w)

    uinit = GridFunction(X)
    for i in range(N):
        uinit.components[i].Set(uin[i], BND, definedon=mesh.Boundaries("bottom"))

    mat_A = CoefficientFunction(tuple([Dcoeff if i==j else 0 for i in range(N) for j in range(N)]),dims=(N,N))

    uprime = CoefficientFunction(tuple([u[i]-u[i]*u[i] if i==j else -u[j]*u[i] for i in range(N) for j in range(N)]),dims=(N,N))

    mat_B = mat_A*uprime

    f =acoeff*u*(1-u)
    normal = specialcf.normal(D)

    for n in range(steps):
        print("Step",n+1,"/",steps)
        a = BilinearForm(X)
        a += SymbolicBFI(- InnerProduct(u,dtv)
                         + (-1)**(n) * eps*InnerProduct(w,v)
                         - (-1)**(n) * InnerProduct(f,v)
                         )
        for j in range(D-1):
            a += SymbolicBFI((-1)**(n) * InnerProduct(mat_B*dxw[j],dxv[j]) )
            a += SymbolicBFI((-1)**(n) * eps*InnerProduct(dxw[j],dxv[j]))

        bottom = 'top' if n%2 else 'bottom'
        top = 'bottom' if n%2 else 'top'
        a += SymbolicBFI((-1)**(n) * InnerProduct(u,v), BND, definedon=mesh.Boundaries(top))
        a += SymbolicBFI((-1)**(n+1) * InnerProduct(CoefficientFunction(tuple(uinit.components[i] for i in range(N))),v), BND, definedon=mesh.Boundaries(bottom))
        if len(nbndc):
            a += SymbolicBFI((-1)**(n+1) * mat_A*InnerProduct(CoefficientFunction(tuple(nbndc[i] for i in range(D-1))),normal[0])*v, BND, definedon=mesh.Boundaries(nbnd))

        gfu = GridFunction(X)

        with TaskManager():
            [success,newtonsteps] = Newton(a,gfu,freedofs=gfu.space.FreeDofs(),maxit=100,maxerr=1e-11,inverse="umfpack",dampfactor=1,printing=True)

        for i in range(N):
            uinit.components[i].Set((1/(1+sum(exp(gfu.components[j]) for j in range(N))))*exp(gfu.components[i]))

    return uinit

if __name__ == "__main__":
    order = 4
    eps=1e-15
    steps=1
    acoeff=1
    C=1
    fu = CoefficientFunction(1/(1+C*exp(-5/6*acoeff*y+1/6*sqrt(6*acoeff)*x))**2)
    nbndc = [CoefficientFunction(-sqrt(2/3*acoeff)*C*exp(-5/6*acoeff*y+1/6*sqrt(6*acoeff)*x)/(1+C*exp(-5/6*acoeff*y+1/6*sqrt(6*acoeff)*x))**3)]
    nbnd="left|right"

    scale=1
    tscale=1

    mesh = Mesh (unit_square.GenerateMesh(maxh=0.5, quad_dominated=False))
    errorold=1
    for h in range(1,7):
        usol = SolFisherKPP(mesh,order,fu,eps,steps,nbndc,nbnd,acoeff,1)
        error=sqrt(Integrate((fu-usol.components[0])**2,mesh))
        rate=log(errorold/error)/log(2)
        errorold=error
        print("===== order: ",order,' error: ', error, ' rate: ', rate)
        input()
        mesh.Refine()
