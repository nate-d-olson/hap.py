# coding=utf-8
#
# Copyright (c) 2010-2016 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

"""
Module for confidence interval calculations, particularly useful for
bioinformatics data analysis in hap.py.
"""

from typing import Dict, Tuple, Union

import numpy as np
import scipy.stats as stats

# Cache for previously calculated values
_VALUE_CACHE: Dict[str, Tuple[float, float, float]] = {}


def jeffreysCI(
    x: Union[int, float], n: Union[int, float], alpha: float = 0.05
) -> Tuple[float, float, float]:
    """Calculate Modified Jeffreys confidence interval for binomial proportions.

    Implements the method from:
    Brown, Cai and DasGupta: Interval Estimation for a Binomial Proportion.
    2001, doi:10.1214/ss/1009213286

    Args:
        x: Number of successes
        n: Number of trials
        alpha: Significance level (default: 0.05 for 95% confidence)

    Returns:
        Tuple of (proportion, lower_bound, upper_bound)
    """
    # HAP-240 avoid division by zero
    if n == 0:
        return 0.0, 0.0, 1.0

    # Check cache for previously calculated values
    key = f"{x}_{n}_{alpha}"
    if key in _VALUE_CACHE:
        return _VALUE_CACHE[key]

    p = x / n
    beta = stats.distributions.beta(x + 0.5, n - x + 0.5)

    # lower bound calculation
    if x == n:
        lower = (alpha / 2) ** (1 / n)
    elif x <= 1:
        lower = 0.0
    else:
        lower = beta.ppf(alpha / 2)

    # upper bound calculation
    if x == 0:
        upper = 1 - (alpha / 2) ** (1 / n)
    elif x >= n - 1:
        upper = 1.0
    else:
        upper = beta.isf(alpha / 2)

    # avoid values outside the unit range due to potential numerical inaccuracy
    lower = max(lower, 0.0)
    upper = min(upper, 1.0)

    # Cache results for future calls
    _VALUE_CACHE[key] = (p, lower, upper)

    return (p, lower, upper)


# Vectorized version for array inputs
binomialCI = np.vectorize(jeffreysCI)
