import polars as pl
import numpy as np
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
from cvxpy import Minimize, Variable, Problem, SCS, ECOS, quad_form, diag, sqrt, norm
import warnings
warnings.filterwarnings("ignore")

def minim(i):
    global sh_data, sh_covinv, sh_Q, sh_mu
    mask = sh_data[:,0] != sh_data[i,0]
    XNi = sh_data[i,1:]
    XM  = sh_data[mask,1:]

    temp = (XNi - XM).T
    x = Variable(mask.sum() , name = "x")

    obj = Minimize(norm(((XNi - x @ XM)@sh_Q).T)+ x @ sqrt(diag(quad_form(temp, sh_covinv))))
    constr  = [x >= 0, x <= 1, sum(x) == 1]

    val = Problem(obj, constr).solve(solver =  ECOS,verbose = False)
    out = x.value.clip(min=0).round(5)
    sh_mu[i, mask] = out
    return i, sh_mu[i]

def init_worker(data, covinv, Q, mu):
    global sh_data, sh_covinv, sh_Q, sh_mu
    sh_data = data
    sh_covinv  = covinv
    sh_Q  = Q
    sh_mu  = mu

if __name__ == "__main__":

    social_vars = ["sexo",
                   "edad_años",
                   "quintil_2",
                   "quintil_3",
                   "quintil_4",
                   "quintil_5",
                   'nivel_escolaridad_primaria',
                   'nivel_escolaridad_secundaria',
                   'nivel_escolaridad_técnica',
                   'nivel_escolaridad_universitaria',
                   "estudia",
                   "trabaja",
                   'macrozona_norte',
                   'macrozona_centro',
                   'macrozona_sur',
                   "horas_trabajo",
                   "n_menores",
                   "n_personas"]

    data = pl.read_csv("data/raw/ENUT_PRE_WEEKEND_IMPUTATION.csv",
                       infer_schema_length=100000,
                       null_values = "NA")
    data =  (
        data
        .with_columns(quintil=pl.col("quintil").cast(int))
        .to_dummies(columns=["quintil", "macrozona", "nivel_escolaridad"])
        .sort("id_persona")
        .select(["dia_fin_semana"] + social_vars ))

    data = data.to_numpy()
    covar = np.cov(data[:,1:].T)

    covinv = np.linalg.inv(covar)
    Q = np.linalg.cholesky(covinv)

    n = len(data)
    mu = np.zeros((n,n))
    with Pool(1, initializer = init_worker, initargs = (data, covinv, Q, mu, )) as p:
        r = list(tqdm(p.imap(minim, range(n)), total=n))

    for i, vec in r:
        mu[i] = vec

    np.save("data/raw/matriz_gemelos2.npy", mu.round(2))
    np.save("data/raw/matriz_gemelos4.npy", mu.round(4))
    np.save("data/raw/matriz_gemelos5.npy", mu.round(5))

    df = pd.DataFrame(mu)
    df.to_csv(
        "data/raw/matriz_gemelos5.csv.gzip",
        header = False,
        index = False,
        compression = "gzip"
    )

    # for i in tqdm(range(n)):
    #     minim(i)
    # np.savetxt("data tesis/matriz_gemelos.txt", mu.round(5), fmt='%.5f',  delimiter=',')
