import numpy as np
from typing import Tuple

def Householder(x: np.array, w: np.array) -> np.array:
	return x - 2 * np.dot(w, x) / np.dot(w, w) * w

def tridiagonalization(A: np.array) -> Tuple[np.array, np.array, np.array]:
	A = A.copy()
	alphas = []
	betas = []

	H = np.identity(np.size(A, 0))

	for m in reversed(range(2, np.size(A, 0))):
		w_i = A[1 :, 0]

		alphas.append(A[0, 0])
		betas.append(np.linalg.norm(w_i))

		w_i[0] -= betas[-1]
		w_i2 = np.dot(w_i, w_i)

		A = A[1 :, 1 :]

		for i in range(m):
			A[:, i] -= 2 * np.dot(w_i, A[:, i]) / w_i2 * w_i

		for row in A:
			row -= 2 * np.dot(w_i, row) / w_i2 * w_i

		for row in H[:, - m :]:
			row -= 2 * np.dot(w_i, row) / w_i2 * w_i

	alphas.extend(np.diag(A))
	betas.append(A[1, 0])

	return (np.array(alphas), np.array(betas), H)

if __name__ == "__main__":
	A = np.array([[2.0, -1.0, 1.0, 3.0], [-1.0, 1.0, 4.0, 2.0], [1.0, 4.0, 2.0, -1.0], [3.0, 2.0, -1.0, 1.0]])
	alphas, betas, H = tridiagonalization(A)
	print(alphas)
	print(betas)
	print(H)