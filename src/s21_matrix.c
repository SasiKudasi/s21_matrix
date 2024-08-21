#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = 0;
  if (result && (rows > 0 && columns > 0)) {
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix) {
      result->rows = rows;
      result->columns = columns;
      for (int i = 0; i < rows && !status; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
        if (!result->matrix[i]) {
          status = 1;
          break;
        }
      }
    } else
      status = 1;
  } else
    status = 1;

  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A && A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i]) free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = 1;
  if (A && B && A->rows == B->rows && A->columns == B->columns && A->matrix) {
    for (int i = 0; i < A->rows && status; i++)
      for (int j = 0; j < A->columns && status; j++)
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07) status = 0;

  } else
    status = 0;
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = 0;
  if (!A || !result || !A->matrix) {
    status = 1;
  } else {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && !status; i++)
      for (int j = 0; j < A->columns && !status; j++) {
        result->matrix[i][j] = number * A->matrix[i][j];
      }
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = 0;
  if (!A || !B || !result || !A->matrix || !B->matrix)
    status = 1;
  else if ((A->rows != B->rows) || (A->columns != B->columns))
    status = 2;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && !status; i++)
      for (int j = 0; j < A->columns && !status; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = 0;
  if (!A || !B || !result)
    status = 1;
  else if (!A->matrix || !B->matrix)
    status = 1;
  else if ((A->rows != B->rows) || (A->columns != B->columns))
    status = 2;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && !status; i++)
      for (int j = 0; j < A->columns && !status; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = 0;
  if (!A || !B || !result)
    status = 1;
  else if (!A->matrix || !B->matrix)
    status = 1;
  else if (A->columns != B->rows)
    status = 2;
  else {
    status = s21_create_matrix(A->rows, B->columns, result) != 0;
    for (int i = 0; i < A->rows && !status; i++)
      for (int j = 0; j < B->columns && !status; j++) {
        for (int k = 0; k < A->columns && !status; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
  }

  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = 0;
  if (!A || !result)
    status = 1;
  else if (!A->matrix)
    status = 1;
  else {
    status = s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows && !status; i++)
      for (int j = 0; j < A->columns && !status; j++)
        result->matrix[j][i] = A->matrix[i][j];
  }
  return status;
}

void get_minor(double **A, double **minor, int skip_row, int ski_column,
               int size) {
  for (int row = 0, x = 0; row < size; row += 1) {
    for (int col = 0, y = 0; col < size; col += 1) {
      if (row != skip_row && col != ski_column) {
        minor[x][y++] = A[row][col];
        if (y == size - 1) {
          y = 0;
          x++;
        }
      }
    }
  }
}

double get_determinant(matrix_t *A, int size) {
  if (size == 1) return A->matrix[0][0];

  matrix_t minor = {0};
  double result = 0;

  s21_create_matrix(size, size, &minor);
  for (int sign = 1, i = 0; i < size; i++, sign *= (-1)) {
    get_minor(A->matrix, minor.matrix, 0, i, size);
    result += sign * A->matrix[0][i] * get_determinant(&minor, size - 1);
  }

  s21_remove_matrix(&minor);
  return result;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = 0;
  if (!A || !result)
    status = 1;
  else if (!A->matrix)
    status = 1;
  else if (A->rows != A->columns)
    status = 2;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    matrix_t minor = {0};

    for (int i = 0; i < A->rows && !status; i++) {
      for (int j = 0; j < A->columns && !status; j++) {
        s21_create_matrix(A->rows - 1, A->rows - 1, &minor);
        get_minor(A->matrix, minor.matrix, i, j, A->rows);
        result->matrix[i][j] =
            ((i + j) % 2 == 0 ? 1 : -1) * get_determinant(&minor, A->rows - 1);
        s21_remove_matrix(&minor);
      }
    }
  }
  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = 0;
  if (!A || !result)
    status = 1;
  else if (!A->matrix)
    status = 1;
  else if (A->rows != A->columns)
    status = 2;
  else
    *result = get_determinant(A, A->rows);

  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det = 0;
  int status = 0;
  if (!A || !result)
    status = 1;
  else if (A->matrix == NULL)
    status = 1;
  else if (A->rows != A->columns)
    status = 2;
  else {
    s21_determinant(A, &det);
    if (det == 0) status = 2;
  }
  if (A->rows == 1 && !status) {
    s21_create_matrix(1, 1, result);
    result->matrix[0][0] = 1.0 / A->matrix[0][0];
  } else if (!status) {
    matrix_t complements = {0};
    matrix_t adjugate = {0};
    s21_calc_complements(A, &complements);
    s21_transpose(&complements, &adjugate);
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && !status; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = adjugate.matrix[i][j] / det;
      }
    }
    s21_remove_matrix(&complements);
    s21_remove_matrix(&adjugate);
  }
  return status;
}