import torch
import numpy as np
from tqdm import tqdm


def find_optimal_PMI_time(M, pi, C, t0=1.0, lr=0.1, epochs=100):
    """
    Find the optimal Markov time for PMI peak with gradient descent.

    :param M: Transition matrix (n * n).
    :param pi: Stationary distribution vector (n).
    :param C: Community membership matrix (n * k).
    :param t0: Initial guess of best Markov time.
    :param lr: Learning rate.
    :param epochs: Number of iterations for optimization.
    :return: the best Markov time.
    """
    M = torch.FloatTensor(M)
    I = torch.eye(M.shape[0])
    Pi = torch.FloatTensor(np.diag(pi))
    pi = torch.FloatTensor(pi)
    C = torch.FloatTensor(C)
    t = torch.nn.Parameter(torch.FloatTensor([t0]))

    optimizer = torch.optim.Adam([t], lr=lr)
    epoch_iterator = tqdm(range(1, epochs+1), desc='Epoch')

    for epoch in epoch_iterator:

        optimizer.zero_grad()
        pmi_matrix = torch.matmul(
            Pi, torch.matrix_exp(-(I-M) * t)
        ) - torch.outer(pi, pi)
        pmi = torch.trace(
            torch.matmul(
                C.t(), torch.matmul(pmi_matrix, C)
            )
        )
        loss = - pmi
        loss.backward()
        optimizer.step()
        t.data.clamp_(min=0)
        epoch_iterator.set_postfix_str(f"t={t.item():.4e}, PMI: {pmi:.4f}")

    return t.item()

