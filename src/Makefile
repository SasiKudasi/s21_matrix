CC = gcc
CC_FLAGS = -Wall -Werror -Wextra -std=c11
FILE = s21_matrix.c
LDFLAGS = -fprofile-arcs -ftest-coverage -lcheck -lm -pthread  -lsubunit

all: s21_matrix.a test 

s21_matrix.a:
	$(CC) $(CC_FLAGS) -c $(FILE)
	ar r s21_matrix.a *.o
	rm *.o
	ranlib s21_matrix.a

test: s21_matrix.a
	$(CC) $(CC_FLAGS) matrix_tests.c s21_matrix.a $(LDFLAGS)
	./a.out
	make clean


gcov: 
	$(CC) $(CC_FLAGS) matrix_tests.c $(FILE) $(LDFLAGS)
	./a.out

gcov_report: gcov
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory coverage_report
	open coverage_report/index.html
	make clean

style:
	clang-format -i --style=Google *.c *.h

clean: 
	rm -rf *.o *.gcno *.out  *.gcda *.info *.a

clean_report:
	rm -rf coverage_report/


rebuild: clean all
