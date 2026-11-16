# -*- coding: utf-8 -*-
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QLabel,
    QComboBox,
    QPushButton,
    QFileDialog,
    QFormLayout,
    QScrollArea,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
import sys
import os

FLEX_SINGLE_PIPETTES = ["flex_1channel_50"]
FLEX_MULTI_PIPETTES = ["flex_8channel_1000"]
FLEX_DEFAULTS = {
    "single_pipette": "flex_1channel_50",
    "multi_pipette": "flex_8channel_1000",
    "single_mount": "right",
    "multi_mount": "left",
    "thermocycler": "thermocyclerModuleV2",
    "magdeck": "magneticBlockV1",
}
FLEX_LABWARE = {
    "24_tuberack_1500ul": ["e14151500starlab_24_tuberack_1500ul"],
    "tiprack_20ul": ["opentrons_flex_96_tiprack_20ul"],
    "tiprack_300ul": ["opentrons_flex_96_tiprack_200ul"],
    "flex_96_tiprack_200ul": ["opentrons_flex_96_tiprack_200ul"],
    "flex_96_tiprack_1000ul": ["opentrons_flex_96_tiprack_1000ul"],
    "clip_source_plate": ["4ti0960rig_96_wellplate_200ul"],
    "clip_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mix_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mag_plate": ["4ti0960rig_96_wellplate_200ul"],
    "final_assembly_plate": ["4ti0960rig_96_wellplate_200ul"],
    "transform_plate": ["4ti0960rig_96_wellplate_200ul"],
    "agar_plate": ["4ti0960rig_96_wellplate_200ul"],
    "12_reservoir_21000ul": ["4ti0131_12_reservoir_21000ul"],
    "96_deepwellplate_2ml": ["4ti0136_96_wellplate_2200ul"],
    "12_corning_wellplate": ["corning_12_wellplate_6.9ml_flat"],
}
ALL_LABWARE_KEYS = list(FLEX_LABWARE)
FLEX_LABWARE_KEYS = set(FLEX_LABWARE.keys())


def unique_options(options):
    seen = []
    for option in options:
        if option not in seen:
            seen.append(option)
    return seen


MODERN_STYLESHEET = "\nQWidget {\n    background-color: #E8F4F8;\n    color: #0A3D62;\n    font-family: 'Helvetica Neue', Arial, sans-serif;\n    font-size: 13px;\n}\n\nQLabel {\n    color: #0A3D62;\n    padding: 2px;\n    background: transparent;\n}\n\n\n\nQComboBox {\n    background-color: #FFFFFF;\n    border: 2px solid #5DADE2;\n    border-radius: 6px;\n    padding: 8px 12px;\n    color: #0A3D62;\n    min-height: 25px;\n    min-width: 200px;\n    font-weight: 500;\n}\n\nQComboBox:hover {\n    border: 2px solid #3498DB;\n    background-color: #F8FBFD;\n}\n\nQComboBox:focus {\n    border: 2px solid #2E86C1;\n}\n\nQComboBox::drop-down {\n    border: none;\n    width: 30px;\n}\n\nQComboBox::down-arrow {\n    image: none;\n    border-left: 5px solid transparent;\n    border-right: 5px solid transparent;\n    border-top: 6px solid #3498DB;\n    margin-right: 8px;\n}\n\n\nQComboBox QAbstractItemView {\n    background-color: #FFFFFF;\n    border: 2px solid #5DADE2;\n    selection-background-color: #AED6F1;\n    selection-color: #0A3D62;\n    padding: 4px;\n}\n\nQPushButton {\n    background-color: #3498DB;\n    color: #FFFFFF;\n    border: none;\n    border-radius: 8px;\n    padding: 12px 24px;\n    font-weight: 600;\n    font-size: 14px;\n    min-height: 20px;\n}\n\nQPushButton:hover {\n    background-color: #2E86C1;\n}\n\nQPushButton:pressed {\n    background-color: #21618C;\n}\n\nQPushButton#generateButton {\n    background-color: #1B4F72;\n    font-size: 16px;\n    padding: 16px 32px;\n    font-weight: bold;\n    letter-spacing: 1px;\n}\n\nQPushButton#generateButton:hover {\n    background-color: #154360;\n}\n\nQScrollArea {\n    border: none;\n    background-color: transparent;\n}\n\nQFrame#headerFrame {\n    background-color: #FFFFFF;\n    border-bottom: 3px solid #3498DB;\n}\n\nQFrame#sectionFrame {\n    background-color: #FFFFFF;\n    border: 2px solid #AED6F1;\n    border-radius: 10px;\n    margin: 8px 0px;\n    padding: 12px;\n}\n\nQLabel#titleLabel {\n    color: #0A3D62;\n    font-size: 28px;\n    font-weight: bold;\n    letter-spacing: 1px;\n}\n\nQLabel#subtitleLabel {\n    color: #5DADE2;\n    font-size: 20px;\n    font-weight: 500;\n}\n\nQLabel#sectionLabel {\n    color: #1B4F72;\n    font-size: 15px;\n    font-weight: bold;\n    padding: 8px 0px;\n}\n\n\nQLabel#fileLabel {\n    color: #3498DB;\n    background-color: #F4F9FD;\n    border: 1px solid #AED6F1;\n    border-radius: 6px;\n    padding: 10px;\n    font-size: 12px;\n}\n"


class GUI(QWidget):
    def __init__(self, root, user_settings):
        self._app = QApplication.instance()
        if self._app is None:
            self._app = QApplication(sys.argv)
        super().__init__()
        self.root = root
        self.user_settings = user_settings
        self.quit_status = True
        self.setStyleSheet(MODERN_STYLESHEET)
        self.setWindowTitle("DNABOT Flex Configuration")
        self.setMinimumWidth(750)
        self.setMinimumHeight(700)
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        container = QWidget()
        self.layout = QVBoxLayout(container)
        self.layout.setContentsMargins(20, 20, 20, 20)
        self.layout.setSpacing(15)
        header = self.create_header()
        self.layout.addWidget(header)
        self.layout.addSpacing(10)
        scroll.setWidget(container)
        main_layout.addWidget(scroll)
        self.add_section_header("Robot Configuration")
        self.layout.addWidget(QLabel("DNABOT Flex"))
        self.add_section_header("Hardware Setup")
        hardware_frame = self.create_frame()
        hardware_layout = QVBoxLayout(hardware_frame)
        self.single_pipette = self.add_dropdown_to_layout(
            hardware_layout, "Single Pipette", FLEX_SINGLE_PIPETTES
        )
        self.multi_pipette = self.add_dropdown_to_layout(
            hardware_layout, "Multi Pipette", FLEX_MULTI_PIPETTES
        )
        self.single_mount = self.add_dropdown_to_layout(
            hardware_layout, "Single Pipette Mount", ["left", "right"]
        )
        self.multi_mount = self.add_dropdown_to_layout(
            hardware_layout, "Multi Pipette Mount", ["left", "right"]
        )
        self.thermocycler = self.add_dropdown_to_layout(
            hardware_layout, "Thermocycler", ["thermocyclerModuleV2"]
        )
        self.magdeck = self.add_dropdown_to_layout(
            hardware_layout, "Magnetic Module", ["magneticBlockV1"]
        )
        self.layout.addWidget(hardware_frame)
        self.add_section_header("Clip Parameters")
        params_frame = self.create_frame()
        params_layout = QVBoxLayout(params_frame)
        params = self.user_settings["parameters"]
        self.premix_linkers = self.add_dropdown_to_layout(
            params_layout, "Premix Linkers", ["Yes", "No"]
        )
        self.premix_linkers.setCurrentText(str(params["premix_linkers"]["value"]))
        self.premix_parts = self.add_dropdown_to_layout(
            params_layout, "Premix Parts", ["Yes", "No"]
        )
        self.premix_parts.setCurrentText(str(params["premix_parts"]["value"]))
        self.linkers_volume = self.add_line_edit_to_layout(
            params_layout, "Linkers Volume (uL)", str(params["linkers_volume"]["value"])
        )
        self.parts_volume = self.add_line_edit_to_layout(
            params_layout, "Parts Volume (uL)", str(params["parts_volume"]["value"])
        )
        self.thermo_temp = self.add_line_edit_to_layout(
            params_layout, "Thermocycler Temp (C)", str(params["thermo_temp"]["value"])
        )
        self.layout.addWidget(params_frame)
        self.add_section_header("Labware Configuration")
        labware_frame = self.create_frame()
        labware_layout = QVBoxLayout(labware_frame)
        self.labware_widgets = {}
        for key in ALL_LABWARE_KEYS:
            items = self.get_labware_options(key)
            widget = self.add_dropdown_to_layout(
                labware_layout, key.replace("_", " ").title(), items
            )
            configured_choice = (
                self.user_settings.get("labwares", {}).get(key, {}).get("id", "")
            )
            if configured_choice and configured_choice in items:
                widget.setCurrentText(configured_choice)
            self.labware_widgets[key] = widget
        self.layout.addWidget(labware_frame)
        self.add_section_header("File Selection")
        file_frame = self.create_frame()
        file_layout = QVBoxLayout(file_frame)
        self.construct_path = ""
        self.source_paths = []
        btn_construct = QPushButton("Select Construct CSV File")
        btn_construct.clicked.connect(self.select_construct)
        file_layout.addWidget(btn_construct)
        self.construct_label = QLabel("No construct file selected")
        self.construct_label.setObjectName("fileLabel")
        self.construct_label.setWordWrap(True)
        file_layout.addWidget(self.construct_label)
        btn_sources = QPushButton("Select BASIC Part CSV Files")
        btn_sources.clicked.connect(self.select_sources)
        file_layout.addWidget(btn_sources)
        self.sources_label = QLabel("No BASIC part files selected")
        self.sources_label.setObjectName("fileLabel")
        self.sources_label.setWordWrap(True)
        self.sources_label.setMinimumHeight(50)
        self.sources_label.setMaximumHeight(150)
        self.sources_scroll = QScrollArea()
        self.sources_scroll.setWidgetResizable(True)
        self.sources_scroll.setWidget(self.sources_label)
        self.sources_scroll.setMaximumHeight(150)
        file_layout.addWidget(self.sources_scroll)
        self.layout.addWidget(file_frame)
        self.layout.addSpacing(20)
        generate_btn = QPushButton("GENERATE PROTOCOL")
        generate_btn.setObjectName("generateButton")
        generate_btn.clicked.connect(self.generate)
        self.layout.addWidget(generate_btn)
        generate_btn.setMaximumWidth(400)
        generate_btn.setMinimumWidth(300)
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(generate_btn)
        button_layout.addStretch()
        self.layout.addLayout(button_layout)
        self.layout.addStretch()
        self.apply_hardware_settings()
        self.show()
        self._app.exec()

    def get_setting_options(self, section, key):
        setting = self.user_settings.get(section, {}).get(key, {})
        options = setting.get("options", [])
        if not isinstance(options, list):
            options = [options] if options else []
        current_value = setting.get("id")
        if current_value and current_value not in options:
            options = [current_value] + options
        return unique_options(
            [str(option) for option in options if option not in (None, "")]
        )

    def get_labware_options(self, key):
        """Return every option declared by this app's settings file."""
        yaml_options = self.get_setting_options("labwares", key)
        if yaml_options:
            return yaml_options
        return unique_options(FLEX_LABWARE.get(key, []))

    def refresh_labware_options(self):
        visible_keys = FLEX_LABWARE_KEYS
        for (key, widget) in self.labware_widgets.items():
            is_visible = key in visible_keys
            widget.setVisible(is_visible)
            widget._label_widget.setVisible(is_visible)
            if not is_visible:
                continue
            current_choice = widget.currentText()
            configured_choice = (
                self.user_settings.get("labwares", {}).get(key, {}).get("id", "")
            )
            items = self.get_labware_options(key)
            widget.clear()
            widget.addItems(items)
            if current_choice and current_choice in items:
                widget.setCurrentText(current_choice)
            elif configured_choice and configured_choice in items:
                widget.setCurrentText(configured_choice)
            elif items:
                widget.setCurrentIndex(0)

    def create_header(self):
        """Create the header with logo and title"""
        header_frame = QFrame()
        header_frame.setObjectName("headerFrame")
        header_frame.setFixedHeight(200)
        header_layout = QVBoxLayout(header_frame)
        header_layout.setContentsMargins(30, 20, 30, 20)
        header_layout.setSpacing(10)
        logo_label = QLabel()
        logo_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "images", "dnabot_logo.png"
        )
        pixmap = QPixmap(logo_path)
        scaled_pixmap = pixmap.scaled(
            200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation
        )
        logo_label.setPixmap(scaled_pixmap)
        logo_label.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(logo_label)
        title = QLabel("Automated BASIC DNA Assembly Protocols")
        title.setObjectName("titleLabel")
        title.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(title)
        subtitle = QLabel("DNABOT 2.3 · Flex")
        subtitle.setObjectName("subtitleLabel")
        subtitle.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(subtitle)
        return header_frame

    def create_frame(self):
        """Create a styled frame for sections"""
        frame = QFrame()
        frame.setObjectName("sectionFrame")
        return frame

    def add_section_header(self, text):
        """Add a section header label"""
        label = QLabel(text)
        label.setObjectName("sectionLabel")
        self.layout.addWidget(label)

    def add_dropdown(self, label, options, callback=None):
        """Add a dropdown to the main layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")
        widget = QComboBox()
        widget.addItems(options)
        if callback:
            widget.currentTextChanged.connect(callback)
        form.addRow(label_widget, widget)
        self.layout.addLayout(form)
        return widget

    def add_dropdown_to_layout(self, parent_layout, label, options, callback=None):
        """Add a dropdown to a specific layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")
        widget = QComboBox()
        widget.addItems(options)
        if callback:
            widget.currentTextChanged.connect(callback)
        widget._label_widget = label_widget
        form.addRow(label_widget, widget)
        parent_layout.addLayout(form)
        return widget

    def add_line_edit_to_layout(self, parent_layout, label, value):
        """Add a line edit to a specific layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")
        widget = QLineEdit()
        widget.setText(value)
        form.addRow(label_widget, widget)
        parent_layout.addLayout(form)
        return widget

    def apply_hardware_settings(self):
        """Update hardware defaults and refresh visible labware rows for the selected robot."""
        self.single_pipette.clear()
        self.single_pipette.addItems(FLEX_SINGLE_PIPETTES)
        self.single_pipette.setCurrentIndex(0)
        self.multi_pipette.clear()
        self.multi_pipette.addItems(FLEX_MULTI_PIPETTES)
        self.multi_pipette.setCurrentIndex(0)
        self.single_pipette.setCurrentText(
            self.user_settings["hardware"]["single_pipette"]["id"]
        )
        self.multi_pipette.setCurrentText(
            self.user_settings["hardware"]["multi_pipette"]["id"]
        )
        self.single_mount.setCurrentText(
            self.user_settings["hardware"]["single_pipette_mount"]["id"]
        )
        self.multi_mount.setCurrentText(
            self.user_settings["hardware"]["multi_pipette_mount"]["id"]
        )
        self.thermocycler.setCurrentText(
            self.user_settings["hardware"]["thermocycler"]["id"]
        )
        self.magdeck.setCurrentText(self.user_settings["hardware"]["mag_deck"]["id"])
        self.refresh_labware_options()

    def select_construct(self):
        """Open file dialog for construct CSV"""
        (file, _) = QFileDialog.getOpenFileName(self, "Select Construct CSV")
        if file:
            self.construct_path = file
            filename = os.path.basename(file)
            self.construct_label.setText(f"✓ Selected: {filename}\n{file}")
        else:
            self.construct_label.setText("No construct file selected")

    def select_sources(self):
        """Open file dialog for source CSV files"""
        (files, _) = QFileDialog.getOpenFileNames(self, "Select BASIC Part CSV Files")
        if files:
            self.source_paths = files
            file_list = "<b>✓ Selected files:</b><br><ul style='margin-top: 5px;'>"
            for file in files:
                filename = os.path.basename(file)
                file_list += f"<li>{filename}</li>"
            file_list += "</ul>"
            self.sources_label.setText(file_list)
        else:
            self.sources_label.setText("No BASIC part files selected")

    def generate(self):
        if (
            not os.path.isfile(self.construct_path)
            or not self.source_paths
            or not all(os.path.isfile(path) for path in self.source_paths)
        ):
            QMessageBox.warning(
                self,
                "Missing input",
                "Select a construct CSV and source CSV files before generating.",
            )
            return
        try:
            linker_volume = float(self.linkers_volume.text())
            part_volume = float(self.parts_volume.text())
            temperature = float(self.thermo_temp.text())
            if linker_volume < 20 or part_volume < 20 or not 4 <= temperature <= 99:
                raise ValueError()
        except ValueError:
            QMessageBox.warning(
                self,
                "Invalid parameters",
                "Premix volumes must be at least 20 uL. Temperature must be between 4 and 99 C.",
            )
            return
        if self.single_mount.currentText() == self.multi_mount.currentText():
            QMessageBox.warning(
                self, "Invalid mounts", "Choose different mounts for the two pipettes."
            )
            return
        """Save settings and close"""
        us = self.user_settings
        us["hardware"]["single_pipette"]["id"] = self.single_pipette.currentText()
        us["hardware"]["multi_pipette"]["id"] = self.multi_pipette.currentText()
        us["hardware"]["thermocycler"]["id"] = self.thermocycler.currentText()
        us["hardware"]["mag_deck"]["id"] = self.magdeck.currentText()
        us["hardware"]["single_pipette_mount"]["id"] = self.single_mount.currentText()
        us["hardware"]["multi_pipette_mount"]["id"] = self.multi_mount.currentText()
        us["parameters"]["premix_linkers"]["value"] = self.premix_linkers.currentText()
        us["parameters"]["premix_parts"]["value"] = self.premix_parts.currentText()
        us["parameters"]["linkers_volume"]["value"] = float(self.linkers_volume.text())
        us["parameters"]["parts_volume"]["value"] = float(self.parts_volume.text())
        us["parameters"]["thermo_temp"]["value"] = float(self.thermo_temp.text())
        for (key, widget) in self.labware_widgets.items():
            us["labwares"][key]["id"] = widget.currentText()
        us["construct_path"] = self.construct_path
        us["sources_paths"] = self.source_paths
        self.quit_status = False
        self.close()
