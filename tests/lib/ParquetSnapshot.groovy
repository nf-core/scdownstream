import groovy.json.JsonOutput
import java.security.MessageDigest

class ParquetSnapshot {
    static Map content(String fragment, def table) {
        [file: new File(fragment).name, schema: table.schema, codecs: table.compressionCodecs, rows_sha256: sha256(rows(table))]
    }

    static Map compactContent(String fragment, def table) {
        [file: new File(fragment).name, codecs: table.compressionCodecs, schema_and_rows_sha256: sha256([table.schema, rows(table)])]
    }

    private static List rows(def table) {
        def columns = table.columnNames
        table.rows.collect { row -> columns.collect { column -> row[column] } }
    }

    private static String sha256(def value) {
        MessageDigest.getInstance('SHA-256').digest(JsonOutput.toJson(value).getBytes('UTF-8')).encodeHex().toString()
    }
}
