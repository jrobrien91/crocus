import panel as pn
pn.extension()

template = pn.template.FastListTemplate(
    sidebar=["# Sidebar Column"],
    main=["# Main Column"],
    header=["# Header Column"],
)
template.servable()
