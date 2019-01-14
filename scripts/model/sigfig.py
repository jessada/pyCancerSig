import numpy as np
#import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.ioff()
import seaborn as sns
from matplotlib.colors import ListedColormap
from collections import defaultdict
from collections import OrderedDict
from matplotlib import ticker
from matplotlib.font_manager import FontProperties
fig_font = FontProperties(family='monospace', weight='bold')

VARIANT_TYPE_COL_NAME = "variant type"
VARIANT_INFO_COL_NAME = "variant info"


def get_profile_dict(feature_infos, feature_fractions):
    profile_dict = defaultdict(dict)
    raw_variant_infos = feature_infos[VARIANT_INFO_COL_NAME]
    raw_variant_types = feature_infos[VARIANT_TYPE_COL_NAME]
    # regrouping data
    for mut_idx in range(len(feature_fractions)):
        variant_info = raw_variant_infos[mut_idx]
        variant_type = raw_variant_types[mut_idx]
        variant_fraction = feature_fractions[mut_idx]
        profile_dict[variant_type][variant_info] = variant_fraction
    return profile_dict

def __get_plottable_variant_info_and_counts(profile_dict):
    """
    Return the tumor context counts in the order in which they are to be plotted.
    """

    # transform dict data into list
    variant_types = list(profile_dict.keys())
    variant_types.sort(key=lambda item: (len(item), item))
    if "MSI" in variant_types:
        # resort variant_types by move MSI to the end
        variant_types.remove("MSI")
        variant_types.append("MSI")
    variant_infos = []
    variant_fractions = []
    variant_type_counts = OrderedDict()
    for variant_type in variant_types:
        raw_variant_infos = sorted(profile_dict[variant_type].keys())
        for variant_info in raw_variant_infos:
            variant_fractions.append(profile_dict[variant_type][variant_info])
        variant_infos += raw_variant_infos
        variant_type_counts[variant_type] = len(raw_variant_infos)
    return variant_type_counts, variant_infos, variant_fractions

def __plot_distribution(variant_type_counts,
                        variant_infos,
                        variant_fractions,
                        title,
                        y_max=None,
                        ):
    """Plot substitution fraction per mutation context"""
    # Set up several subplots
    n_variant_types = len(variant_type_counts)
    fig, axes = plt.subplots(nrows=1, ncols=n_variant_types, figsize=(20, 5.5))
    fig.canvas.set_window_title(title)
    fig.suptitle(title, fontsize=14)

    if y_max is None:
        y_max = max(variant_fractions) * 1.2
#    y_max = 0.05
    y_min = min(variant_fractions) * 1.2

    # Set up some colors and markers to cycle through...
    colors = sns.color_palette("bright", n_colors=len(variant_type_counts)).as_hex()

    graph = 0
    first_variant_idx = 0
    for ax, variant_count, color in zip(axes, variant_type_counts.values(), colors):
        x = np.arange(variant_count) - 10
        fractions = variant_fractions[first_variant_idx:first_variant_idx+variant_count]
        ax.bar(x, fractions, color=color, width=.5, align='edge')
        # Labels for the rectangles
        new_ticks = variant_infos[first_variant_idx:first_variant_idx+variant_count]
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, (end - start) / (float(variant_count)*1.1)) + .85)
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(new_ticks))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, font_properties=fig_font, color='k')
        # Standardize y-axis ranges across subplots (in percentage units)
        ax.set_ylim([y_min, y_max])
        vals = ax.get_yticks()
        ax.set_yticklabels(['{:3.0f}%'.format(val * 100) for val in vals])
        plt.setp(ax.yaxis.get_majorticklabels(), color='k', fontweight='bold')
        graph += 1
        first_variant_idx += variant_count

    # Set labels
    axes[0].set_ylabel('Mutation Type Probability', fontweight='bold', color='k')

    labels = variant_type_counts.keys()

    for ax, label in zip(axes, labels):
        ax.set_xlabel(label, fontweight='bold', color='k', size=13)

    # Remove boundaries between subplots
    for ax in axes[1:]:
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.set_yticklabels([])
    axes[0].spines['right'].set_color('none')

    # Leave off last x-tick to reduce clutter.
    for ax in axes:
        xticks = ax.get_xticks()
        ax.set_xticks(xticks[0:-1])

    # Merge subplots together so that they look like one graph
    fig.subplots_adjust(wspace=-.03)
    return y_max

def plot_distribution(title,
                      profile_dict=None,
                      feature_infos=None,
                      feature_fractions=None,
                      y_max=None,
                      ):
    if profile_dict is None:
        profile_dict = get_profile_dict(feature_infos, feature_fractions)
    variant_type_counts, variant_infos, variant_fractions = __get_plottable_variant_info_and_counts(profile_dict)
    return __plot_distribution(variant_type_counts,
                               variant_infos,
                               variant_fractions,
                               title,
                               y_max=y_max,
                               )
