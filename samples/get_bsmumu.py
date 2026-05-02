import ROOT as r 
import pandas as pd

def get_bsmumu():
    files = [ 
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000010_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000009_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000001_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000003_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000008_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000004_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000007_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000005_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000006_1.tuple.root",
        "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2025/TUPLE.ROOT/00376532/0000/00376532_00000002_1.tuple.root"
    ]

    df = r.RDataFrame("MCDT/MCDecayTree", files) 
    node = df.Range(500000)
    print("snapshotting mcdt")
    node.Snapshot("MCDecayTree", "bsmumu_rdf.root")
    print("done")
def parquet( ):
    import pandas as pd
    df = r.RDataFrame("MCDecayTree", "bsmumu_rdf.root")
    dfPd = pd.DataFrame(df.AsNumpy())
    dfPd.to_parquet("bsmumu_rdf.parquet")

def plot1d_histo(
    ax: plt.Axes,
    df: pd.DataFrame,
    column: str,
    query: Optional[str] = None,
    weights: Optional[str] = None,
    density: bool = False,
    label: Optional[str] = None,
    **hist_kwargs: Any,
) -> None:
    """
    Plot a 1D histogram on a given Axes from a DataFrame.
    """
    data = df.query(query) if query is not None else df
    
    w = data[weights] if weights is not None else None
    
    ax.hist(
        data[column],
        weights=w,
        density=density,
        label=label,
        **hist_kwargs,
    )
    ax.set_xlabel(column)
    if density:
        ax.set_ylabel("Density")
    else:
        ax.set_ylabel("Entries")


def plot2d_histo(
    ax: plt.Axes,
    df: pd.DataFrame,
    x: str,
    y: str,
    query: Optional[str] = None,
    bins: int | Tuple[int, int] = 50,
    range: Optional[Tuple] = None,
    cmin: float = 1,
    label: Optional[str] = None,
    **hist2d_kwargs: Any,
) -> None:
    """
    Plot a 2D histogram (colormesh) on a given Axes.
    """
    data = df.query(query) if query is not None else df

    h, xedges, yedges = np.histogram2d(
        data[x], data[y], bins=bins, range=range
    )
    
    # Mask empty bins
    h = np.ma.masked_where(h < cmin, h)
    
    im = ax.pcolormesh(xedges, yedges, h.T, **hist2d_kwargs)
    
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    plt.colorbar(im, ax=ax, label="Entries" if label is None else label)


def plot1d_efficiency(
    ax: plt.Axes,
    df: pd.DataFrame,
    num: str,           # numerator column (e.g. "passed")
    den: str,           # denominator column (e.g. "all")
    var: str,           # variable to bin on
    query: Optional[str] = None,
    bins: int | np.ndarray = 30,
    range: Optional[Tuple[float, float]] = None,
    label: Optional[str] = None,
    **step_kwargs: Any,
) -> None:
    """
    Plot efficiency (num/den) vs a variable as a step histogram.
    Useful for turn-on curves, selection efficiencies, etc.
    """
    data = df.query(query) if query is not None else df
    
    # Bin the variable
    hist_num, bin_edges = np.histogram(data[var][data[num]], bins=bins, range=range)
    hist_den, _        = np.histogram(data[var][data[den]],  bins=bins, range=range)
    
    # Avoid division by zero
    eff = np.divide(
        hist_num, hist_den, 
        out=np.zeros_like(hist_num, dtype=float), 
        where=hist_den > 0
    )
    
    # Plot as step function
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.step(bin_centers, eff, where='mid', label=label, **step_kwargs)
    
    ax.set_xlabel(var)
    ax.set_ylabel("Efficiency")
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)    
if __name__ == "__main__":
    # get_bsmumu()
    
    parquet() 