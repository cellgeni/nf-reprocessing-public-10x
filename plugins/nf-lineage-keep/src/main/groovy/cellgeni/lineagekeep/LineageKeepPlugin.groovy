package cellgeni.lineagekeep

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Entry point of the nf-lineage-keep plugin.
 */
@CompileStatic
class LineageKeepPlugin extends BasePlugin {

    LineageKeepPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }
}
