/*
 * Copyright 2015 Jaroslaw Glowacki
 *
 * This file is part of cgDNArecon.
 *
 * cgDNArecon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cgDNArecon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cgDNArecon.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * The testing application that runs all the tests.
 */

#include<cppunit/TestFixture.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <algebra3d/global.h>
#include <iostream>

#include "Tester.h"

using namespace std;

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

	bool wasSuccessful = runner.run("", false);
	return !wasSuccessful;
}
