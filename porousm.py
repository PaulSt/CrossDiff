import sys
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import SegMesh
from netgen.geom2d import unit_square
from ngsolve.solvers import *
from ngsolve import *

def SolPorousM(mesh,order,uin,m,eps,steps,nbndc=[],nbnd="",shift=0,acoeff=1):
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

    mat_A = CoefficientFunction(tuple([acoeff*m*pow(u[i],m-1) if i==j else 0 for i in range(N) for j in range(N)]),dims=(N,N))

    uprime = CoefficientFunction(tuple([u[i]-u[i]*u[i] if i==j else -u[j]*u[i] for i in range(N) for j in range(N)]),dims=(N,N))

    mat_B = mat_A*uprime

    normal = specialcf.normal(D)

    for n in range(shift,steps):
        print("Step",n+1,"/",steps)
        a = BilinearForm(X)
        a += SymbolicBFI(- (-1)**(n) * InnerProduct(u,dtv)
                         + eps*InnerProduct(w,v)
                         )
        for j in range(D-1):
            a += SymbolicBFI(InnerProduct(mat_B*dxw[j],dxv[j]) )
            a += SymbolicBFI(eps*InnerProduct(dxw[j],dxv[j]))

        A=5
        B=5
        nbndc=[CoefficientFunction((x-A)/(m*(m+1)*(B-(y+n))))]
        bottom = 'top' if n%2 else 'bottom';
        top = 'bottom' if n%2 else 'top';
        a += SymbolicBFI(InnerProduct(u,v), BND, definedon=mesh.Boundaries(top))
        a += SymbolicBFI(-InnerProduct(uin,v), BND, definedon=mesh.Boundaries(bottom))
        a += SymbolicBFI(-mat_A*nbndc[0]*normal[0]*v, BND, definedon=mesh.Boundaries("left||right"))

        gfu = GridFunction(X)

        with TaskManager():
            Newton(a,gfu,freedofs=gfu.space.FreeDofs(),maxit=100,maxerr=1e-11,inverse="umfpack",dampfactor=1,printing=True)

        for i in range(N):
            uinit.components[i].Set((1/(1+sum(exp(gfu.components[j]) for j in range(N)))) * (exp(gfu.components[i])))

    return uinit


if __name__ == "__main__":
    m=2
    eps=0
    steps=1
    acoeff=1
    k=(m-1)/acoeff
    fu = CoefficientFunction(pow(k*x+k*y+0.5,1/(m-1)))
    A=5
    B=5
    fu = CoefficientFunction(((x-A)*(x-A)/(2*m*(m+1)*(B-y))))
    nbndc=[CoefficientFunction((x-A)/(m*(m+1)*(B-y)))]
    nbnd="left||right"
    shift=0
    scale=1
    tscale=1
    order=4
    errorold=1
    for h in range(1,7):
        mesh = Mesh (unit_square.GenerateMesh(maxh=1/2**h, quad_dominated=False))

        usol = SolPorousM(mesh,order,fu,m,eps,steps,nbndc,nbnd,0,acoeff)
        error=sqrt(Integrate((fu-usol.components[0])**2,mesh))
        rate=log(errorold/error)/log(2)
        errorold = error
        print("===== order: ",order,' error: ', error, ' rate: ', rate)
        input()
