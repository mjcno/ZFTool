**ZFtool** is a toolbox for **quantification of cellular proliferation in vivo** of the **zebrafish** trough image processing methods.

Zebrafishes are injected with a *green fluorescence protein* (*GFP*) in order to visualize the cancer mass.

ZFTool eliminates the auto-fluorescence of the zebrafish through computation of area with different intensity thresholds and automatically computing the auto-fluorescence threshold, which is established for both images at 0hpi and 24-48-72hpi (_hpi=hours post injection_), respectively. ZFTool then computes the cancer mass area and mean intensity for both images and calculates the **proliferation index**.

Zebrafish (_Danio rerio_) is a model organism that has emerged as a tool for **cancer
research**, cancer being the second most common cause of death after cardiovascular
disease for humans in the developed world.

Zebrafish is a useful model for
xenotransplantation of human cancer cells and toxicity studies of different
chemotherapeutic compounds in vivo. Compared to the murine model, the zebrafish
model is **faster**, can be screened using **high-throughput methods** and has a **lower
maintenance cost**, making it possible and affordable to create personalized therapies.

While several methods for cell proliferation determination based on image acquisition
and quantification have been developed, some drawbacks still remain. In the
xenotransplantation technique, quantification of cellular proliferation in vivo is critical to
standardize the process for future preclinical applications of the model.

**ZFtool can establish a base threshold that eliminates embryo auto-fluorescence and
measures the area of marked cells (GFP) and the intensity of those cells to define a
*proliferation index*.**

## Usage

The analyzed images must be presented in **their own folder** following this format, being `N` the fish number and `XX` the hours elapsed after injection (hpi):

 * Grey image at 0hpi, named `N - 0h bn.tif`
 * Grey image at 24, 48 or 72 hpi, named `N - XXh bn.tif`
 * GFP fluorescence image at 0hpi, named `N - 0h gfp.tif`
 * GFP fluorescence image at 24, 48 or 72 hpi, named `N - XXh bn.tif`

There is an example included in the repository, in the folder `Images`:

* `2 - 0h bn.tif`
* `2 - 0h gfp.tif`
* `2 - 48h bn.tif`
* `2 - 48h gfp.tif`


With the included images you can run the following example:

 * Image folder: `Images`
 * Fish number: `2`
 * Initial measurement hour (0hpi): `0`
 * Final measurement hour (24hpi, 48hpi or 72hpi): `48`
