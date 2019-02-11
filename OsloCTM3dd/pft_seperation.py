# Coniferous forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_cf>=95)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax41, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_cf>=30)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax41, label=labels[j], color=colors[j], ls=linestyle)
# Decidous forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_df>=70)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax42, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_df>=54)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax42, label=labels[j], color=colors[j], ls=linestyle)
# Tropical forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_tf>=98)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax43, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_tf>=98)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax43, label=labels[j], color=colors[j], ls=linestyle)
# Cropland
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_ac>=70)&(dd_velo[j].lat>=40),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax44, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_ac>=50)&(dd_velo[j].lat<=-30),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax44, label=labels[j], color=colors[j], ls=linestyle)
# Grassland
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_gr>=63)&(dd_velo[j].lat>=40),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax45, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_gr>=74)&(dd_velo[j].lat<=-30),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax45, label=labels[j], color=colors[j], ls=linestyle)
# Tundra
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_tu>=75)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax46, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_tu>=22)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax46, label=labels[j], color=colors[j], ls=linestyle)
# Snow and ice (some other barren ground?)
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_is>=98)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax47, label=labels[j], color=colors[j], ls=linestyle)
            print(j, (dd_velo[j].where((pft_is>=98)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).mean())
        else:
            (dd_velo[j].where((pft_is>=98)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax47, label=labels[j], color=colors[j], ls=linestyle)
# Ocean
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_oc>=98)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax48, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_oc>=98)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax48, label=labels[j], color=colors[j], ls=linestyle)
# Desert
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        if b_nh:
            (dd_velo[j].where((pft_de>=98)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax49, label=labels[j], color=colors[j], ls=linestyle)
        else:
            (dd_velo[j].where((pft_de>=93)&(dd_velo[j].lat<0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax49, label=labels[j], color=colors[j], ls=linestyle)
