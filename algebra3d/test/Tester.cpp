/*
 * Copyright 2015 Jaroslaw Glowacki
 * jarek (dot) glowacki (at) gmail (dot) com
 *
 * This file is part of algebra3d.
 *
 * algebra3d is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * algebra3d is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with algebra3d.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * The testing application that runs all the tests.
 */
#include<cstdlib>
#include<cppunit/TestFixture.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iostream>
#include <algebra3d/global.h>

#include "Tester.h"

using namespace std;

bool DO_HEAVY_TESTS = false;

int main(int argc, char **argv) {
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry =
			CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest(registry.makeTest());

	if (sizeof(fpType) == sizeof(float)) {
		cout << "Testing single precission" << endl;
	} else {
		cout << "Testing double precission" << endl;
	}
	cout << "Precision of comparision: " << PRECISION << endl;

	if(argc > 1 && atoi(argv[1]) > 0) {
		DO_HEAVY_TESTS = true;
	}

	bool wasSuccessful = runner.run("", false);
	return !wasSuccessful;
}
