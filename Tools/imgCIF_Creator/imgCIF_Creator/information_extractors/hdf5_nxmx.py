"""Functionalities to extract information from an hdf5 file.
"""

import sys
from collections import defaultdict
import h5py as h5
import numpy as np
from imgCIF_Creator.command_line_interfaces import parser
from . import extractor_interface, extractor_utils
from imgCIF_Creator.output_creator import imgcif_creator


class Extractor(extractor_interface.ExtractorInterface):
    """See also documentation of the init method.

    Args:
        extractor_interface (class): the interface that must be implemented by
            every extractor class
    """

    def __init__(self, filename) -> None:
        """This extractor allows to extract the scan and setup information from
        hdf5 NeXus NxMx files. When an instance of the extractor is initialized
        then the information is attempted to be extracted and stored in class
        attributes. The extractor provides public methods to make this information
        accessible.

        Args:
            filename (str): the filename of the hdf5 master file
        """

        super().__init__()

        file_mapping = self._get_scan_file_mapping(filename)
        self.scan_info = self._get_scan_info(filename, file_mapping)
        self._first_scan = sorted(self.scan_info.keys())[0]
        n_frames_first_scan = self.scan_info[self._first_scan][1]['frames']
        frames_per_file = self._get_frames_per_file(
            file_mapping, n_frames_first_scan)
        self.all_frames = self._get_all_frames(file_mapping, frames_per_file)

        self.setup_info = self._get_setup_info(filename, file_mapping)
        self.axes = self.get_axes_description(filename, file_mapping)


    def get_misc_info(self):
        """Return the information that was found about the doi and the array
        intensities overload.

        Returns:
            dict: a dictionary containing the array intensities overload
        """

        doi = None
        overload = self.setup_info.get('overload')
        # TODO can there be an hdf5 entry with the temp?
        temperature = self.setup_info.get('temperature')

        return {'doi' : doi,
                'overload' : overload,
                'temperature' : temperature,
                }


    def get_source_info(self):
        """Return the information about the facility and beamline or the instrument,
        model and location. Cif block: _diffrn_source

        Returns:
            dict: a dictionary containing the information about the source
        """

        facility = self.setup_info.get('facility')
        beamline = self.setup_info.get('beamline')
        manufacturer = None
        model = self.setup_info.get('detector_model')
        location = None

        source_info = {'beamline' : beamline,
                       'facility' : facility,
                       'manufacturer' : manufacturer,
                       'model' : model,
                       'location' : location}

        return source_info


    def get_axes_info(self):
        """Return the information about the axes settings. Cif block: _axis

        Returns:
            dict: a dictionary containing the information about the axes settings
        """

        axes, axis_type, equip, depends_on, vector, offset = \
            self._make_axes(self.axes)

        axes_info = {'axes' : axes,
                     'axis_type' : axis_type,
                     'equip' : equip,
                     'depends_on' : depends_on,
                     'vector' : vector,
                     'offset' : offset}

        return axes_info


    def get_array_info(self):
        """Return the information about the array. Cif block: _array_structure_list_axis
        and _array_structure_list

        Returns:
            dict: a dictionary containing the information about the array
        """

        pixel_size = \
            (self.setup_info.get('x_pixel_size'), self.setup_info.get('y_pixel_size'))

        if self.setup_info.get('fast_direction') == "horizontal":
            array_precedence = [1, 2]
        else:
            array_precedence = [2, 1]

        array_info = {
            'axis_id' : None,
            'axis_set_id': None,
            'pixel_size' : pixel_size,
            'array_id' : None,
            'array_index' : None,
            'array_dimension' : self.setup_info.get('array_dimensions'),
            'array_direction' : self.setup_info.get('fast_direction'),
            'array_precedence' : array_precedence,
        }

        return array_info


    def get_detector_info(self):
        """Return the information about the detector. Cif block: _diffrn_detector
        and _diffrn_detector_axis.

        Returns:
            dict: a dictionary containing the information about the detector
        """

        #TODO multiple detectors are not supportet (yet?)
        detector_info = {
            'detector_id' : None,
            'number_of_axes' : None,
            'axis_id' : None,
            'detector_axis_id' : None
        }

        return detector_info


    def get_radiation_info(self):
        """Return the information about the wavelength an type of radiation.
        Cif block: _diffrn_radiation and _diffrn_radiation_wavelength

        Returns:
           dict: a dictionary containing the information about the radiation
        """

        # TODO not creating a list of wavelength id's here (yet)
        wavelength = self.scan_info[self._first_scan][1]['wavelength']
        rad_type = self.scan_info[self._first_scan][1]['rad_type']

        return {'rad_type' : rad_type,
                'wavelength' : wavelength}


    def get_scan_settings_info(self):
        """Return the information about the scans, this is a dictionary containing
        the starting point settings of the axes and the details of each scan.

        For example for scan '08':
        {'08': ({'chi': -60.991, 'phi': 110.0, 'detector_2theta': -12.4,
        'omega': -18.679, 'distance': 40.0}, {'frames': 12, 'axis': 'omega',
        'incr': 2.0, 'time': 1800.0, 'start': -40.679, 'range': 24.0,
        'wavelength': 0.560834, 'x_pixel_size': 0.172, 'y_pixel_size': 0.172,
        'mini_header': ['# detector: pilatus100k',.... ],
        ...}

        Returns:
            dict: a dictionary containing the information about the scans
        """

        return self.scan_info


    def _get_scan_info(self, master_file, file_mapping):
        """Retrieve the scan information from the master file.

        Args:
            master_file (str): the filename of the master file
            file_mapping (dict): a mapping between the scan and the files (if this
                are multiple files the values are lists)

        Raises:
            Exception: scan range does not match increment

        Returns:
            dict: a dictionary containing the scan information
        """

        scan_info = {}
        goniometer_rot_direction =  \
            parser.CommandLineParser().request_input('goniometer_rot_direction')

        with h5.File(master_file, 'r') as h5_master:
            for scan_no, scan in file_mapping.keys():
                scan_axis_found = False
                start_axes_settings = {}
                h5item = h5_master.get(f'{scan}/sample')

                sample = h5_master.get(f'{scan}/sample')
                for path in [grp.name for grp in sample.values()]:
                    h5item = h5_master.get(path)
                    if not isinstance(h5item, h5.Group) or len(h5item.keys()) > 1:
                        continue
                    axis = list(h5item.keys())[0]
                    axis = self._replace_names(axis)
                    # drop sam axes
                    if axis in ['sam_x', 'sam_y', 'sam_z']:
                        continue
                    # take only first
                    dataset = list(h5item.values())[0]
                    try:
                        if len(dataset) > 1:
                            print(f'\nIdentified {axis} as scan axis')
                            scan_axis_found = True
                            scan_axis = axis
                            n_frames = len(dataset)
                            scan_start = dataset[0]
                            scan_stop = dataset[-1]
                            path = f'{scan}/sample/transformations/{scan_axis}_increment_set'
                            scan_incr = self._get_hdf5_item(h5_master, path)[0]
                            # if the rotation direction is the other way round,
                            # change increment
                            if goniometer_rot_direction in ['anticlockwise', 'a']:
                                scan_incr = scan_incr * -1
                                scan_stop *= -1
                            # because of 0.1*137 = 13.700000000000001 we round
                            scan_range = round(scan_incr * n_frames, 10)
                    except TypeError:
                        # for fast and slow this is a scalar and no list
                        pass

                    # TODO convert units if neccessary?
                    start_axes_settings[axis] = dataset[0]

                path = f'{scan}/instrument/detector_z/det_z'
                h5item = self._get_hdf5_item(h5_master, path)
                start_axes_settings['trans'] = h5item[0]

                if not scan_axis_found:
                    print('No scan axis found! Aborting!')
                    sys.exit()

                # Check increment and range match
                if not np.isclose(scan_start + scan_incr * (n_frames - 1), scan_stop, atol=1e-6):
                    raise Exception(
                        f"Scan range does not match increment: \
    {scan_start} to {scan_stop}, {n_frames} steps of {scan_incr}")

                path = f'{scan}/instrument/detector/count_time'
                exposure = self._get_hdf5_item(h5_master, path)
                path = f'{scan}/instrument/beam/incident_wavelength'
                wavelength = self._get_hdf5_item(h5_master, path)
                path = f'{scan}/instrument/source/type'
                rad_type = self._get_hdf5_item(h5_master, path)
                scan_details = {"frames" : n_frames,
                                "axis" : scan_axis,
                                "incr" : scan_incr,
                                "time" : exposure,
                                "start" : scan_start,
                                "range" : scan_range,
                                "wavelength" : wavelength,
                                "rad_type" : rad_type,
                                }
                scan_info[scan_no] = (start_axes_settings, scan_details)

        scan_info = extractor_utils.prune_scan_info(scan_info)

        return scan_info


    def _get_all_frames(self, file_mapping, frames_per_file):
        """Return a dictionary containing all frames mapped with their corresponding
        filename and hdf5 path.

        Args:
            file_mapping (dict): a mapping between the scan and the files (if this
                are multiple files the values are lists)
            frames_per_file (list): a list containing the number of frames per file

        Returns:
            dict: a dictionary with the scan number and the total frame number as
                key (tuple) and a dictionary with filename, hdf5 path and frame
                in file as keys.
        """

        all_frames = {}
        for scan_no, scan in sorted(file_mapping):
            counter = 1
            for i, entry in enumerate(sorted(file_mapping[(scan_no, scan)])):
                filename, dpath = entry
                for frame_in_file in range(1, frames_per_file[i] + 1):
                    all_frames[(scan_no, counter)] = {
                        'filename' : filename,
                        'path' : dpath,
                        'frame' : frame_in_file
                    }
                    counter += 1

        return all_frames


    def _get_scan_file_mapping(self, file_name):
        """Obtain the frames from the master.

        Args:
            file_name (str): the file name

        Returns:
            defaultdict(list): a dictionary with a list of the external file links
        """

        scan_file_mapping = defaultdict(list)
        with h5.File(file_name, 'r') as h5_master:
            # if there are multiple scans the ordering is not ensured
            for scan_no, scan in enumerate(h5_master.keys(), 1):
                group = h5_master[f'{scan}/data']
                for key in group:
                    link = group.get(key, getlink=True)
                    if isinstance(link, h5.ExternalLink):
                        scan_file_mapping[(scan_no, scan)].append(
                            (link.filename, link.path))

        return scan_file_mapping


    def _get_frames_per_file(self, file_mapping, n_frames):
        """Return the frames per file from the user input and ensure that the number
        of frames agrees with the total number of frames.

        Args:
            file_mapping (dict): a mapping between the scan and the files (if this
                are multiple files the values are lists)
            n_frames (int): the total number of frames

        Returns:
            list: a list with the number of frames per file
        """

        first_key = list(file_mapping.keys())[0]
        n_files = len(file_mapping[first_key])
        frame_numbers = []
        while len(frame_numbers) != n_files:
            # strangely this in encapsulated into two lists
            files = [tup[0] for tup in list(file_mapping.values())[0]]
            print(f"\nFound {n_files} files with {n_frames} frames in total: \n\
 {', '.join(files)}\n", end='')
            frame_numbers = parser.CommandLineParser().request_input('frame_numbers')
            frame_numbers = frame_numbers.replace(' ', '').split(',')
            frame_numbers = [int(number) for number in frame_numbers]
            if sum(frame_numbers) != n_frames:
                print(f'The sum of the frames per file is not matching the sum in \
the master file ({sum(frame_numbers)} != {n_frames}), please try again.')
                frame_numbers = []
            elif len(frame_numbers) != n_files:
                print('The numbers of frames per file does not match the number \
of files. Please try again!')

        return frame_numbers


    def _get_setup_info(self, master_file, file_mapping):
        """Return the setup infomration that is found in the hdf5 file.

        Args:
            master_file (str): the filename of the hdf5 master file
            file_mapping (dict): a mapping between the scan and the files (if this
                are multiple files the values are lists)

        Returns:
            dict: a dictionary containing the information retrieved from hdf5
        """

        setup_info = {}
        with h5.File(master_file, 'r') as h5_master:
            for _, scan in file_mapping.keys():

                path = f'{scan}/instrument/source/name'
                setup_info['facility'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/source/beamline' #can this path in exist?
                setup_info['beamline'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/entry/instrument/detector/description'
                setup_info['detector_model'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/detector/saturation_value'
                setup_info['overload'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/detector/x_pixel_size'
                setup_info['x_pixel_size'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/detector/y_pixel_size'
                setup_info['y_pixel_size'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/detector/module/data_size'
                setup_info['array_dimensions'] = self._get_hdf5_item(h5_master, path)

                path = f'{scan}/instrument/detector/module/fast_pixel_direction'
                fast_direction = self._get_hdf5_item(h5_master, path, 'vector')

                if fast_direction[0] != 0 and fast_direction[1] == fast_direction[2] == 0:
                    setup_info['fast_direction'] = 'horizontal'
                elif fast_direction[1] != 0 and fast_direction[0] == fast_direction[2] == 0:
                    setup_info['fast_direction'] = 'vertical'

        return setup_info


    def get_axes_description(self, master_file, file_mapping):
        """Retrieve the goniometer axis settings from hdf5.

        Args:
            master_file (str): the filename of the master file
            file_mapping (dict): a mapping between the scan and the files (if this
                are multiple files the values are lists)

        Returns:
            axes (dict): a dictionary containing the gonimeter axes and
                their vectors, offsets, etc.
        """

        axes = {}
        with h5.File(master_file, 'r') as h5_master:
            for axis in imgcif_creator.GONIOMETER_AXES + imgcif_creator.DETECTOR_AXES:
                axis = axis.lower()
                scan = list(file_mapping.keys())[0][1]
                if axis in imgcif_creator.GONIOMETER_AXES:
                    path = f'{scan}/sample/sample_{axis}/{axis}'
                elif axis in imgcif_creator.DETECTOR_AXES:
                    if axis == 'det_z':
                        path = f'{scan}/instrument/detector_z/{axis}'
                    else:
                        path = f'{scan}/instrument/detector/module/{axis}'

                h5item = h5_master.get(path)
                if h5item is None or len(h5item.attrs) == 0:
                    continue

                axis_type = h5item.attrs['transformation_type'].decode('utf-8')
                depends_on = h5item.attrs['depends_on'].decode('utf-8').split('/')[-1]
                # loop through a dependency chain until it does not depend
                # on sam_x etc. anymore
                depends_on_h5item = h5_master[h5item.attrs['depends_on']]
                while depends_on in ['sam_x', 'sam_y', 'sam_z', 'module_offset']:
                    depends_on = depends_on_h5item.attrs['depends_on']\
                        .decode('utf-8').split('/')[-1]
                    depends_on_h5item = h5_master[depends_on_h5item.attrs['depends_on']]
                    depends_on_h5item = h5_master[list(depends_on_h5item.attrs.items())\
                        [0][1].decode('utf-8')]

                depends_on = self._replace_names(depends_on)
                axis = self._replace_names(axis)
                if axis == 'dety':
                    depends_on = 'detx'
                elif axis == 'trans':
                    depends_on = 'two_theta'

                vector = h5item.attrs['vector']
                axes[axis] = \
                    {
                        'depends_on': depends_on,
                        'vector': vector,
                        'axis_type': axis_type
                    }

                # TODO is an omega axis always present?
                if axis == 'omega':
                    # two_theta axis is duplicate of omega except for the equipment
                    axes['two_theta'] = \
                        {
                            'depends_on': depends_on,
                            'vector': vector,
                            'axis_type': axis_type
                        }

            # conversion from NeXus/McStas to CBF coordinate system
            # the kind of transformations (rotation) which have to be performed depend
            # on the position and rotation direction of the goniometer (seen from
            # the source)
            # TODO is this correct?
            if axes['omega']['vector'].all() == np.array([-1, 0, 0]).all():
                goniometer_pos = 'right'
            elif axes['omega']['vector'].all() == np.array([1, 0, 0]).all():
                goniometer_pos = 'left'
            else:
                print('Could not identify goniometer position! Aborting!')
                sys.exit()
            print(f'\nIdentified goiniometer position (from source): {goniometer_pos}')

            for axis in axes:
                axes[axis]['offset'] = np.array([0.0, 0.0, 0.0])

            path = f'{scan}/instrument/detector/module/module_offset'
            axes['detx']['offset'] = self._get_hdf5_item(h5_master, path, 'offset')

            path = f'{scan}/instrument/detector_distance'
            # off_z = self._get_hdf5_item(h5_master, path)
            # off_z = off_z[0] if isinstance(off_z, np.ndarray) else off_z
            axes['trans']['offset'] = np.array([0.0, 0.0, 0.0])

        for axis, content in axes.items():
            axes[axis]['vector'] = \
                self._rotate_from_nexus_to_cbf(content['vector'], goniometer_pos)
            axes[axis]['offset'] = \
                self._rotate_from_nexus_to_cbf(content['offset'], goniometer_pos)

        return axes


    def _replace_names(self, name):
        """Replace some axis names with the correct one.

        Args:
            name (str): the axis name to replace

        Returns:
            str: the axis name replaced
        """

        # this is hardcoded
        if name == 'fast_pixel_direction':
            return 'detx'
        if name == 'slow_pixel_direction':
            return 'dety'
        if name == 'det_z':
            return 'trans'

        return name


    def _rotate_from_nexus_to_cbf(self, vector, goniometer_pos):
        """Rotate the vectors from NeXus to CBF coordinate system.

        Args:
            vector (string): The vector that shall be rotated as string e.g. '1 0 0'
            goniometer_pos (string): The position of the goniometer seen from the
                source ('left' or 'right')

        Returns:
            rotated_vector (string): The rotated vector.
        """

        if goniometer_pos == 'right':
            rotated_vector = self._rotate_vector(vector, 'y')
        elif goniometer_pos == 'left':
            rotated_vector = self._rotate_vector(vector, 'x')

        return rotated_vector


    def _rotate_vector(self, vector, rotation_axis):
        """Rotate a vector around a rotation axis.

        Args:
            vector (string): The vector that shall be rotated as string e.g. '1 0 0'
            rotation_axis (string): The axis around which shall be rotated
                'x', 'y' or 'z')

        Returns:
            rotated_vector (string): The rotated vector.
        """

        if rotation_axis == 'x':
            rot = vector * np.array([1, -1, -1])
        elif rotation_axis == 'y':
            rot = vector * np.array([-1, 1, -1])
        elif rotation_axis == 'z':
            rot = vector * np.array([-1, -1, 1])
        else:
            rot = None

        return rot

    def _get_hdf5_item(self, h5_master, path, attr=None):
        """Return the hdf5 item at the given path with optional attribute. Return
        the utf-8 decoded string if possible, return None if not present.

        Args:
            h5_master (h5py._hl.files.File): the hdf5 master file opened with h5py
            path (str): the path where the item should be found
            attr (str, optional): the optional attribute of the hdf5 entry. Defaults to None.

        Returns:
            any: the value obtained from the hdf5 file (string, float, list, etc.)
        """

        h5item = h5_master.get(path)
        if h5item is not None:
            if attr is not None:
                h5item = h5item.attrs[attr]
            h5item = h5item[()]
            try:
                h5item = h5item.decode('utf-8')
            except AttributeError:
                pass

        return h5item


    def _make_axes(self, axes_info):
        """Return the axes information in the required format for get_axes_info.

        Args:
            axes (dict): a dictionary containing the gonimeter axes and
                their vectors, offsets, etc.

        Returns:
            list: a list with the axis properties in the format for get_axes_info
        """

        axes = []
        axis_type = []
        equip = []
        depends_on = []
        vector = []
        offset = []

        for axis in axes_info:
            axes.append(axis)
            ax_type = axes_info[axis]['axis_type']
            axis_type.append(ax_type)
            if axis in imgcif_creator.GONIOMETER_AXES:
                equip.append('goniometer')
            elif axis in imgcif_creator.DETECTOR_AXES:
                equip.append('detector')
            else:
                equip.append('equipment')
            depends_on.append(axes_info[axis]['depends_on'])
            vector.append(list(axes_info[axis]['vector']))
            offset.append(list(axes_info[axis]['offset']))

        return [axes, axis_type, equip, depends_on, vector, offset]
