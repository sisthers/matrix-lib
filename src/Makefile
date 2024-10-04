CC = g++ -Werror -Wextra -Wall -std=c++17
HEADERS = *.h
SRCS = console_ui.cc
TESTS = matrix_tests.cc
TESTEX = matrix_tests
PKG = `pkg-config --cflags --libs gtest`

all: gcov_report
	
gcov_report: clean
	$(CC) -fprofile-arcs -ftest-coverage $(SRCS) $(TESTS) -o gcovreport \
	$(PKG)
	./gcovreport
	lcov -t "gcovreport" --no-external -o gcovreport.info -c -d .
	genhtml -o report gcovreport.info
	open report/index.html

test: clean $(LIB)
	$(CC) $(TESTS) $(LIB) -o $(TESTEX) \
	$(PKG)
	./$(TESTEX)

fix_style:
	clang-format -i $(SRCS) $(TESTS) $(HEADERS)

style:
	clang-format -n $(SRCS) $(TESTS) $(HEADERS)

leaks: clean $(LIB)
	$(CC) $(TESTS) $(LIB) -o $(TESTEX) \
	$(PKG)
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(TESTEX)

check: style leaks clean

clean:
	rm -rf *.o *.a *.gcno *.gcda *.info $(TESTEX) gcovreport report 

rebuild: clean all
