from decimal import getcontext as decimalcontext

import pytest


@pytest.fixture(params=[5, 10, 100, 1000])
def precision(request):
    prec = request.param
    decimalcontext().prec = prec + 1
    return prec
