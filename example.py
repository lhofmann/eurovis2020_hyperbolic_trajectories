#!/usr/bin/env python
#
# Command line arguments:
#   input vector field (e.g., vector_field.vti)
#   name of the point data velocity array (e.g., u)
#   input candidate line (e.g., candidate_line.vtp)
#
# If no arguments are given, the skewing oscillating gyre saddle and a analytical candidate line are used instead.
# The discretized input data is written to vector_field.vti and candidate_line.vtp.
#

from __future__ import division

from vtk import (
    vtkDataObject,
    vtkImageData,
    vtkLineSource,
    vtkPolyData,
    vtkXMLImageDataReader,
    vtkXMLImageDataWriter,
    vtkXMLPolyDataReader,
    vtkXMLPolyDataWriter,
)
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import numpy as np
import sys
from hyperbolic_trajectories import vtkApproximateDHT


def main():
    dht = vtkApproximateDHT()

    if len(sys.argv) == 4:
        vector_field_reader = vtkXMLImageDataReader()
        vector_field_reader.SetFileName(sys.argv[1])
        candidate_line_reader = vtkXMLPolyDataReader()
        candidate_line_reader.SetFileName(sys.argv[3])

        dht.SetInputConnection(0, vector_field_reader.GetOutputPort())
        dht.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, sys.argv[2])
        dht.SetInputConnection(1, candidate_line_reader.GetOutputPort())
    else:
        # construct analytical vector field
        vector_field = vtkSkewingOscillatingGyreSaddle()
        # output sampled vector field
        writer = vtkXMLImageDataWriter()
        writer.SetInputConnection(vector_field.GetOutputPort(0))
        writer.SetFileName('vector_field.vti')
        writer.Update()
        # construct analytical candidate line
        line_source = vtkLineSource()
        line_source.SetResolution(100)
        line_source.SetPoint1(0., 0., 0.1)
        line_source.SetPoint2(0., 0., 7.9)
        line_source.Update()
        candidate_line = line_source.GetOutput()
        points = vtk_to_numpy(candidate_line.GetPoints().GetData())
        t = points[:, 2]
        points[:, 0] = 0.2 * np.sin(t * (2.0 * np.pi) / 4.0)
        points[:, 1] = 0.8 * np.sin(t * (2.0 * np.pi) / 4.0)
        # output sampled candidate line
        writer = vtkXMLPolyDataWriter()
        writer.SetInputDataObject(candidate_line)
        writer.SetFileName('candidate_line.vtp')
        writer.Update()

        dht.SetInputConnection(0, vector_field.GetOutputPort(0))
        dht.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, 'u')
        dht.SetInputDataObject(1, candidate_line)

    writer = vtkXMLPolyDataWriter()
    writer.SetInputConnection(dht.GetOutputPort())
    writer.SetFileName('dht.vtp')
    writer.Update()


class vtkSkewingOscillatingGyreSaddle(VTKPythonAlgorithmBase):
    def __init__(self):
        self._dimensions = [100, 100, 100]
        self._bounds = [-2.5, 2.5, -2.5, 2.5, 0., 8.]
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkImageData')

    def RequestInformation(self, request, inInfo, outInfo):
        executive = self.GetExecutive()
        out_info = outInfo.GetInformationObject(0)
        extent = [0, self._dimensions[0] - 1, 0, self._dimensions[1] - 1, 0, self._dimensions[2] - 1]
        out_info.Set(executive.WHOLE_EXTENT(), extent, 6)
        return 1

    def velocity(self, x, y, t):
        def gyre_saddle(a, x, y):
            def clamped_cos(x):
                return np.cos(np.clip(x, -np.pi / 2.0, np.pi / 2.0))

            def inner_saddle(x, y):
                dx = -np.sin(np.pi * x) * np.cos(np.pi * y) + a * np.sin(np.pi * y) * np.cos(np.pi * x)
                dy = np.sin(np.pi * y) * np.cos(np.pi * x) - a * np.sin(np.pi * x) * np.cos(np.pi * y)
                return dx, dy

            def outer_saddle(x, y):
                k = (y >= np.abs(x)).astype(float) - (y <= -np.abs(x)).astype(float)
                l = (x > np.abs(y)).astype(float) - (x < -np.abs(y)).astype(float)
                dx = k * a * clamped_cos(np.pi * x - a * np.pi * (y - k / 2)) \
                     - l * clamped_cos(np.pi * y - a * np.pi * (x - l / 2))
                dy = k * clamped_cos(np.pi * x - a * np.pi * (y - k / 2)) \
                     - l * a * clamped_cos(np.pi * y - a * np.pi * (x - l / 2))
                return dx, dy

            inside = (np.abs(x) <= 0.5) & (np.abs(y) <= 0.5)
            dx_inner, dy_inner = inner_saddle(x, y)
            dx_outer, dy_outer = outer_saddle(x, y)
            dx = np.where(inside, dx_inner, dx_outer)
            dy = np.where(inside, dy_inner, dy_outer)
            return dx, dy

        a = .5 * (np.sin(t * (2. * np.pi) / 1.0))

        x -= 0.2 * np.sin(t * (2.0 * np.pi) / 4.0)
        y -= 0.8 * np.sin(t * (2.0 * np.pi) / 4.0)

        dx, dy = gyre_saddle(a, x, y)
        dx *= 1.5
        dy *= 1.5
        dt = np.ones_like(t)

        return np.stack([dx, dy, dt], axis=-1)

    def RequestData(self, request, inInfo, outInfo):
        executive = self.GetExecutive()
        out_info = outInfo.GetInformationObject(0)
        image_output = vtkImageData.GetData(outInfo, 0)

        coords = list(np.linspace(self._bounds[2 * i], self._bounds[2 * i + 1], self._dimensions[i]) 
                      for i in range(3))
        grid = np.meshgrid(*coords, indexing='ij')

        u = self.velocity(*grid)

        image_output.SetDimensions(self._dimensions)
        image_output.SetOrigin(self._bounds[::2])
        image_output.SetSpacing([coords[i][1] - coords[i][0] for i in range(3)])

        u_cont = np.ascontiguousarray(u.reshape((-1, u.shape[-1]), order='F'))
        array = numpy_to_vtk(u_cont, deep=True)
        array.SetName('u')
        image_output.GetPointData().AddArray(array)
        image_output.GetPointData().SetActiveVectors('u')

        return 1


if __name__ == '__main__':
    main()
