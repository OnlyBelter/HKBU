import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error


def get_interaction_matrix(r1, r2, a12, a21, k1, k2):
    return np.array([[r1/k1, r1*a12/k1],
                     [r2*a21/k2, r2/k2]])


def dx_dt(x, t=0, r=None, a=None, k=None):
    """
    Return the growth rate of fox and rabbit populations.
    # competitive LV model, not predator-prey model
    # https://en.wikipedia.org/wiki/Competitive_Lotka%E2%80%93Volterra_equations
    """
    x[x < 1e-9] = 0
    k = k.reshape(-1, 1)
    x = x.reshape(-1, 1)
    r = r.reshape(-1, 1)
    # print(f'>>> {t}: {np.ravel(x * (r + int_max @ x))}')
    # print(x.shape, a.shape, k.shape)
    return np.ravel(x * r * (1 - a @ x * 1/k))


def forward_predict(dx_dt, x_0, t, r, a, k):
    x_pred, infodict = integrate.odeint(dx_dt, x_0, t, args=(r, a, k), full_output=True)
    print(infodict['message'])
    return x_pred


def obj_func(ra, t, k, x_true, dx_dt):
    r, a = to_matrix(ra)
    x_0 = x_true[0, :]
    x_pred = integrate.odeint(dx_dt, x_0, t, args=(r, a, k))
    return mean_squared_error(y_true=x_true, y_pred=x_pred)


def dX_dt2(x, t=0, r=None):
    # int_max = np.array([[10, 7, 12],
    #                     [15, 10, 8],
    #                     [7, 11, 10]])
    int_max = np.array([[10, 7],
                        [15, 10]])
    return np.ravel(x * (r - int_max @ x))


def generate_data(n_species=3, k=None, t=None, result_file_path=None):
    n = n_species
    # r = np.array([1, 1])  # r1, r2
    a = np.array([[1, 0.7, 1.2],
                  [1.5, 1, 0.8],
                  [0.7, 1.1, 1]])  # alpha12, alpha21
    # t = np.linspace(0, 100, m)  # time
    x_0 = np.random.random(n)  # initialize conditions: n species
    # print(x_0)
    r = np.array([10] * n)  # inherent growth rate
    x = forward_predict(dx_dt=dx_dt, a=a, r=r, x_0=x_0, k=k, t=t)
    # print(len(x), x)
    if result_file_path is not None:
        x_df = pd.DataFrame(data=x, index=t, columns=[f'species{i}' for i in range(n_species)])
        x_df.to_csv(result_file_path, float_format='%g')
    return x


def to_vector(r, a, n_species=3):
    assert r.shape == (n_species, 1)
    assert a.shape == (n_species, n_species)
    return np.hstack([r.flatten(), a.flatten()])


def to_matrix(vec, n_species=3):
    """

    :param vec: r (3 by 1) + a (3 by 3)
    :param n_species:
    :return:
    """
    assert vec.shape == (n_species + n_species ** 2,)
    return vec[:n_species].reshape(n_species, 1), vec[n_species:].reshape(n_species, n_species)


def plot_density_over_time(t, x, result_file_path=None, title=None):
    """

    :param t: m time points
    :param x: m by n_species
    :param result_file_path:
    :return:
    """
    plt.figure(figsize=(8, 6))
    for i, species in enumerate(x.T):
        # species1, species2, species3 = X.T
        plt.plot(t, species, label=f'species{i+1}')
        # p.plot(t, species2, 'b-', label='species2')
        # p.plot(t, species3, 'y-', label='species3')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('population')
    if title is not None:
        plt.title(title)
    if result_file_path is not None:
        plt.savefig(result_file_path, dpi=200)


if __name__ == '__main__':
    # Definition of parameters
    # np.random.seed(42)  # fix random seed for reproducibility
    n_species = 3
    m = 1000
    t = np.linspace(0, 100, m)  # time
    k = np.array([1] * n_species)  # k1, k2， carrying capacity, fix to 1 for all species
    result_fp = 'raw_data.csv'
    if not os.path.exists(result_fp):
        # x_true, m by n_species
        x_true = generate_data(result_file_path=result_fp, t=t, k=k, n_species=n_species)
    else:
        print(f'Using generated dataset: {result_fp}')
        x_true = pd.read_csv(result_fp, index_col=0)
        x_true = x_true.values
    plot_density_over_time(t=t, x=x_true, result_file_path='3species_competitive_LV_model_x_true.png',
                           title='Ground truth')

    # estimate parameters of a and r
    # x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    r = np.random.random((n_species, 1))
    a = np.random.random((n_species, n_species))
    res = minimize(obj_func, to_vector(r=r, a=a), args=(t, k, x_true, dx_dt), tol=1e-8)
    r_pred, a_pred = to_matrix(res.x)
    print(to_matrix(res.x))
    x_pred = forward_predict(dx_dt=dx_dt, a=a_pred, r=r_pred, x_0=x_true[0, :], k=k, t=t)
    plot_density_over_time(t=t, x=x_pred, result_file_path='3species_competitive_LV_model_x_pred.png',
                           title='Predicted r and a')
